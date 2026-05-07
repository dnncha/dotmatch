import type { Metadata } from "next";
import "./globals.css";

export const metadata: Metadata = {
  title: "DotMatch - Exact Short-DNA Assignment",
  description:
    "DotMatch is a fast exact short-DNA known-target assignment engine for CRISPR guides, barcodes, primers, panels, and whitelists.",
  openGraph: {
    title: "DotMatch",
    description:
      "Exact one-edit known-target assignment with deterministic ambiguity semantics and workflow-ready FASTQ outputs.",
    type: "website"
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
