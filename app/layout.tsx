import type { Metadata } from "next";
import "./globals.css";

const siteUrl = process.env.NEXT_PUBLIC_SITE_URL ?? "https://dnncha.github.io/dotmatch";
const socialImageUrl = `${siteUrl}/dotmatch-og.png`;
const twitterImageUrl = `${siteUrl}/dotmatch-twitter.png`;

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
    type: "website",
    url: siteUrl,
    images: [
      {
        url: socialImageUrl,
        width: 1200,
        height: 630,
        alt: "DotMatch exact known-target short-DNA assignment"
      }
    ]
  },
  twitter: {
    card: "summary_large_image",
    title: "DotMatch",
    description:
      "Exact one-edit known-target assignment with deterministic ambiguity semantics and workflow-ready FASTQ outputs.",
    images: [
      {
        url: twitterImageUrl,
        width: 1200,
        height: 630,
        alt: "DotMatch exact known-target short-DNA assignment"
      }
    ]
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
