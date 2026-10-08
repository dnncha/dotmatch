#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

static int mutated;

static void mutate_destination(void) {
    const char *path = getenv("DOTMATCH_MUTATE_PATH");
    const char *mode = getenv("DOTMATCH_MUTATE_MODE");
    if (mutated || path == NULL || mode == NULL) return;
    mutated = 1;
    const char payload[] = "external bytes\n";
    char temporary[PATH_MAX];
    int replace = strcmp(mode, "replace") == 0;
    int fd;
    if (replace) {
        snprintf(temporary, sizeof(temporary), "%s.concurrent-XXXXXX", path);
        fd = mkstemp(temporary);
    } else {
        fd = open(path, O_WRONLY | (strcmp(mode, "create") == 0 ? O_CREAT | O_EXCL : 0), 0600);
    }
    if (fd < 0) return;
    struct stat original;
    int have_original = fstat(fd, &original) == 0;
    if (write(fd, payload, sizeof(payload) - 1) < 0) { close(fd); return; }
    if (strcmp(mode, "edit") == 0 && have_original) {
        struct timespec times[] = {original.st_atim, original.st_mtim};
        futimens(fd, times);
    }
    close(fd);
    if (replace) {
        int (*real_rename)(const char *, const char *) = dlsym(RTLD_NEXT, "rename");
        real_rename(temporary, path);
    }
}

int fclose(FILE *stream) {
    const char *trigger = getenv("DOTMATCH_MUTATE_TRIGGER");
    if (!mutated && trigger != NULL && strcmp(trigger, "close") == 0) {
        char descriptor[64], path[PATH_MAX];
        snprintf(descriptor, sizeof(descriptor), "/proc/self/fd/%d", fileno(stream));
        ssize_t length = readlink(descriptor, path, sizeof(path) - 1);
        if (length >= 0) {
            path[length] = '\0';
            if (strstr(path, ".dotmatch-stage-") != NULL) mutate_destination();
        }
    }
    int (*real_fclose)(FILE *) = dlsym(RTLD_NEXT, "fclose");
    return real_fclose(stream);
}

int rename(const char *source, const char *destination) {
    const char *fail_path = getenv("DOTMATCH_FAIL_PUBLICATION_PATH");
    if (fail_path != NULL && strcmp(destination, fail_path) == 0 &&
        strstr(source, ".dotmatch-stage-") != NULL) {
        const char *trigger = getenv("DOTMATCH_MUTATE_TRIGGER");
        if (trigger != NULL && strcmp(trigger, "failure") == 0) mutate_destination();
        errno = EIO;
        return -1;
    }
    int (*real_rename)(const char *, const char *) = dlsym(RTLD_NEXT, "rename");
    return real_rename(source, destination);
}
