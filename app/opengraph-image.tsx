import { ImageResponse } from "next/og";

export const alt = "DotMatch exact known-target short-DNA assignment";
export const size = {
  width: 1200,
  height: 630
};
export const contentType = "image/png";

const bases = "ACGTACGTTGCAAGTCGATCGTACCTAGGCTA".split("");

export default function Image() {
  return new ImageResponse(
    (
      <div
        style={{
          width: "100%",
          height: "100%",
          display: "flex",
          flexDirection: "column",
          justifyContent: "space-between",
          padding: "64px 72px",
          background:
            "linear-gradient(135deg, #f8faf8 0%, #ffffff 52%, #e2f7ef 100%)",
          color: "#101513",
          fontFamily: "Inter, Arial, sans-serif"
        }}
      >
        <div style={{ display: "flex", alignItems: "center", gap: 22 }}>
          <div
            style={{
              width: 70,
              height: 70,
              border: "5px solid #101513",
              borderRadius: 14,
              background:
                "linear-gradient(90deg, transparent 44%, rgba(14,124,90,0.34) 44% 56%, transparent 56%), linear-gradient(0deg, transparent 44%, rgba(29,102,209,0.28) 44% 56%, transparent 56%), #ffffff"
            }}
          />
          <div style={{ display: "flex", flexDirection: "column" }}>
            <span style={{ fontSize: 32, fontWeight: 800 }}>DotMatch</span>
            <span style={{ color: "#58645f", fontSize: 24 }}>
              Apache-2.0 open core for real FASTQ workflows
            </span>
          </div>
        </div>

        <div style={{ display: "flex", flexDirection: "column", gap: 24 }}>
          <div
            style={{
              maxWidth: 930,
              fontSize: 72,
              lineHeight: 0.98,
              fontWeight: 860,
              letterSpacing: 0
            }}
          >
            Exact known-target short-DNA assignment.
          </div>
          <div
            style={{
              maxWidth: 900,
              color: "#33423b",
              fontSize: 32,
              lineHeight: 1.28
            }}
          >
            CRISPR guides, inline barcodes, primers, panels, and whitelists
            with deterministic unique, ambiguous, and no-match semantics.
          </div>
        </div>

        <div
          style={{
            display: "flex",
            alignItems: "center",
            justifyContent: "space-between",
            gap: 30
          }}
        >
          <div style={{ display: "flex", gap: 14 }}>
            {["C core", "CLI", "Python", "Benchmarks", "Claim gates"].map((item) => (
              <span
                key={item}
                style={{
                  padding: "12px 18px",
                  border: "1px solid #dfe7e2",
                  borderRadius: 10,
                  background: "#ffffff",
                  color: "#0e7c5a",
                  fontSize: 22,
                  fontWeight: 760
                }}
              >
                {item}
              </span>
            ))}
          </div>
          <div style={{ display: "flex", gap: 6 }}>
            {bases.map((base, index) => (
              <span
                key={`${base}-${index}`}
                style={{
                  width: 18,
                  height: 38,
                  borderRadius: 5,
                  background:
                    base === "A" || base === "T" ? "#23b082" : "#1d66d1",
                  opacity: index % 5 === 0 ? 1 : 0.54
                }}
              />
            ))}
          </div>
        </div>
      </div>
    ),
    size
  );
}
