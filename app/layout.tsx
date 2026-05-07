import type { Metadata } from "next";
import "./globals.css";

const siteUrl = process.env.NEXT_PUBLIC_SITE_URL ?? "https://donncha.github.io/dotmatch";

export const metadata: Metadata = {
  metadataBase: new URL(siteUrl),
  title: "DotMatch - Exact Short-DNA Assignment",
  description:
    "DotMatch is a fast exact short-DNA known-target assignment engine for CRISPR guides, barcodes, primers, panels, and whitelists.",
  keywords: [
    "bioinformatics",
    "computational biology",
    "CRISPR",
    "FASTQ",
    "barcode demultiplexing",
    "edit distance",
    "known-target assignment"
  ],
  openGraph: {
    title: "DotMatch",
    description:
      "Exact one-edit known-target assignment with deterministic ambiguity semantics and workflow-ready FASTQ outputs.",
    type: "website"
  },
  twitter: {
    card: "summary_large_image",
    title: "DotMatch",
    description:
      "Exact one-edit known-target assignment with deterministic ambiguity semantics and workflow-ready FASTQ outputs."
  }
};

export default function RootLayout({
  children
}: Readonly<{
  children: React.ReactNode;
}>) {
  return (
    <html lang="en">
      <body>{children}</body>
    </html>
  );
}
