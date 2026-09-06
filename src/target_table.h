/* Internal target-table reader. Included by qda.c after seq_table helpers.
 * Keep this input contract aligned with python/dotmatch/target_io.py.
 * gzip is detected by zlib; no compressed bytes become literal target rows.
 */
#ifndef DOTMATCH_TARGET_TABLE_H
#define DOTMATCH_TARGET_TABLE_H

#define TARGET_TABLE_MAX_LINE (1024U * 1024U)
#define TARGET_TABLE_MAX_FIELDS 128U

static int target_table_line(gzFile input, char **buf, size_t *capacity) {
    size_t used = 0;
    int c;
    while ((c = gzgetc(input)) != -1) {
        if (c == '\0' || used >= TARGET_TABLE_MAX_LINE) return -1;
        if (used + 2 > *capacity) {
            size_t next = *capacity ? *capacity * 2 : 1024;
            char *grown = (char *)realloc(*buf, next);
            if (grown == NULL) return -1;
            *buf = grown;
            *capacity = next;
        }
        (*buf)[used++] = (char)c;
        if (c == '\n') break;
    }
    if (c == -1) {
        int error = Z_OK;
        (void)gzerror(input, &error);
        if (error != Z_OK && error != Z_STREAM_END) return -1;
    }
    if (!used) return 0;
    while (used && ((*buf)[used - 1] == '\n' || (*buf)[used - 1] == '\r')) --used;
    (*buf)[used] = '\0';
    return 1;
}

static void target_table_trim(char *text) {
    size_t length = strlen(text), start = 0;
    while (start < length && (text[start] == ' ' || text[start] == '\t')) ++start;
    while (length > start && (text[length - 1] == ' ' || text[length - 1] == '\t')) --length;
    memmove(text, text + start, length - start);
    text[length - start] = '\0';
}

/* Decode one CSV/TSV row in place. Reject unclosed/stray quotes and overflow
 * rather than splitting a quoted identifier or silently discarding columns.
 */
static int target_table_fields(char *line, char delimiter, char **fields) {
    char *read = line, *write = line;
    size_t count = 0;
    for (;;) {
        if (count == TARGET_TABLE_MAX_FIELDS) return -1;
        fields[count++] = write;
        while (*read == ' ' && delimiter != ' ') ++read;
        if (*read == '"') {
            ++read;
            for (;;) {
                if (!*read) return -1;
                if (*read == '"') {
                    if (read[1] == '"') { *write++ = '"'; read += 2; }
                    else { ++read; break; }
                } else { *write++ = *read++; }
            }
            while (*read == ' ') ++read;
            if (*read && *read != delimiter) return -1;
        } else {
            while (*read && *read != delimiter) {
                if (*read == '"') return -1;
                *write++ = *read++;
            }
        }
        int more = *read == delimiter;
        if (more) ++read;
        *write++ = '\0';
        target_table_trim(fields[count - 1]);
        if (!more) break;
    }
    return (int)count;
}

static int target_table_alias(const char *name, const char *const *aliases) {
    for (size_t i = 0; aliases[i]; ++i) if (field_eq(name, aliases[i])) return 1;
    return 0;
}

static int target_table_column(char **fields, size_t count, const char *const *aliases) {
    int found = -1;
    for (size_t i = 0; i < count; ++i) {
        if (!target_table_alias(fields[i], aliases)) continue;
        if (found != -1) return -2;
        found = (int)i;
    }
    return found;
}

static int target_table_valid_text(const char *text, int is_sequence) {
    for (const unsigned char *p = (const unsigned char *)text; *p; ++p) {
        if (*p < 32 || *p == 127 || (is_sequence && (*p == 32 || *p > 127))) return 0;
    }
    return 1;
}

static int read_target_table(const char *path, seq_table *table) {
    static const char *const ids[] = {"target_id", "guide_id", "barcode_id", "id", "name", "sgrna", "guide", "sgrnaid", "sgrna_id", NULL};
    static const char *const seqs[] = {"target_seq", "guide_seq", "barcode_seq", "sequence", "seq", "grna.sequence", "bases", "sgrna.sequence", "sgrna_sequence", "sgrna_seq", "guide_sequence", "guidesequence", NULL};
    static const char *const genes[] = {"gene", "gene_id", "gene_symbol", "gene.symbol", "target_gene", NULL};
    gzFile input = gzopen(path, "rb");
    if (input == NULL) return -1;
    (void)gzbuffer(input, 128 * 1024);
    char *line = NULL;
    size_t capacity = 0, line_number = 0, width = 0;
    int id_col = 0, seq_col = 1, gene_col = 2, named = 0, sequence_only = 0;
    int first = 1, status = 0;
    char delimiter = '\t';
    const char *failure = "invalid target table";
    while ((status = target_table_line(input, &line, &capacity)) == 1) {
        ++line_number;
        if (line_number == 1 && strlen(line) >= 3 && memcmp(line, "\xef\xbb\xbf", 3) == 0)
            memmove(line, line + 3, strlen(line + 3) + 1);
        char *first_char = line;
        while (*first_char == ' ' || *first_char == '\t') ++first_char;
        if (!*first_char || *first_char == '#') continue;
        if (first) {
            size_t n = strlen(path);
            if (n >= 3 && field_eq(path + n - 3, ".gz")) n -= 3;
            int csv_path = n >= 4 && strncasecmp(path + n - 4, ".csv", 4) == 0;
            delimiter = csv_path || (strchr(line, ',') && !strchr(line, '\t')) ? ',' : '\t';
        }
        char *fields[TARGET_TABLE_MAX_FIELDS];
        int parsed = target_table_fields(line, delimiter, fields);
        if (parsed < 1) { failure = "malformed quoted row or more than 128 columns"; goto invalid; }
        size_t nf = (size_t)parsed;
        if (first) {
            int id = target_table_column(fields, nf, ids);
            int seq = target_table_column(fields, nf, seqs);
            int gene = target_table_column(fields, nf, genes);
            if (id == -2 || seq == -2 || gene == -2) { failure = "multiple possible ID, sequence or gene columns; provide an unambiguous library"; goto invalid; }
            sequence_only = nf == 1 && seq == 0;
            named = sequence_only || (id >= 0 && seq >= 0);
            width = nf;
            first = 0;
            if (named) {
                for (size_t i = 0; i < nf; ++i) {
                    if (!fields[i][0]) { failure = "empty header column"; goto invalid; }
                    for (size_t j = 0; j < i; ++j) if (field_eq(fields[i], fields[j])) {
                        failure = "duplicate header columns"; goto invalid;
                    }
                }
                id_col = id; seq_col = seq; gene_col = gene;
                continue;
            }
        }
        if (nf != width) { failure = "row has a different number of fields from the header or first row"; goto invalid; }
        const char *id, *seq, *gene = "";
        char generated[48];
        if (sequence_only || (!named && nf == 1)) {
            (void)snprintf(generated, sizeof(generated), "target_%zu", table->count);
            id = generated; seq = fields[0];
        } else {
            if (id_col < 0 || seq_col < 0 || (size_t)id_col >= nf || (size_t)seq_col >= nf) {
                failure = "missing target columns"; goto invalid;
            }
            id = fields[id_col]; seq = fields[seq_col];
            if (gene_col >= 0 && (size_t)gene_col < nf) gene = fields[gene_col];
        }
        if (id[0] == '"' || gene[0] == '"') {
            failure = "target ID and gene may not start with a literal double quote (TSV identity boundary)"; goto invalid;
        }
        if (!id[0] || !seq[0]) { failure = "target ID and sequence must be non-empty"; goto invalid; }
        if (!target_table_valid_text(id, 0) || !target_table_valid_text(gene, 0) || !target_table_valid_text(seq, 1)) {
            failure = "target text contains control characters or sequence whitespace/non-ASCII"; goto invalid;
        }
        if (push_record_gene(table, id, strlen(id), seq, strlen(seq), gene, strlen(gene)) != 0) {
            failure = "could not allocate target library"; goto invalid;
        }
    }
    if (status < 0) { failure = "corrupt input, embedded NUL, or target row exceeds 1 MiB"; goto invalid; }
    free(line);
    if (gzclose(input) != Z_OK) return -1;
    return table->count ? 0 : -1;
invalid:
    fprintf(stderr, "%s:%zu: %s\n", path, line_number, failure);
    free(line);
    (void)gzclose(input);
    return -1;
}
#endif
