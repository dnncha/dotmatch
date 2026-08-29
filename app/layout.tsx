import type { Metadata, Viewport } from "next";
import "./globals.css";

const siteUrl = process.env.NEXT_PUBLIC_SITE_URL ?? "https://dnncha.github.io/dotmatch";
const socialImageUrl = `${siteUrl}/dotmatch-og.png`;
const twitterImageUrl = `${siteUrl}/dotmatch-twitter.png`;
const socialImageAlt =
  "DotMatch preview showing unique, ambiguous, unmatched, and invalid read assignments";

export const metadata: Metadata = {
  metadataBase: new URL(siteUrl),
  applicationName: "DotMatch",
  title: "DotMatch - Clear Read Assignment for Known DNA Targets",
  description:
    "DotMatch compares a chosen part of each sequencing read with known short DNA targets and keeps uncertain or invalid assignments visible.",
  authors: [{ name: "DotMatch maintainers", url: "https://github.com/dnncha/dotmatch" }],
  creator: "DotMatch maintainers",
  publisher: "DotMatch",
  category: "Bioinformatics software",
  alternates: {
    canonical: siteUrl
  },
  keywords: [
    "bioinformatics",
    "computational biology",
    "CRISPR",
    "FASTQ",
    "known-target sequencing",
    "assignment reliability",
    "barcode demultiplexing",
    "barcode panel design",
    "barcode troubleshooting",
    "FASTQ sequence matching",
    "nf-core",
    "MultiQC",
    "Galaxy workflows",
    "Snakemake",
    "core facility sequencing",
    "bioinformatics workflow adoption",
    "sequencing core facility QC"
  ],
  openGraph: {
    title: "DotMatch - Clear Read Assignment",
    description:
      "See which sequencing reads match known DNA targets and which remain ambiguous, unmatched, or invalid.",
    type: "website",
    siteName: "DotMatch",
    locale: "en_US",
    url: siteUrl,
    images: [
      {
        url: socialImageUrl,
        secureUrl: socialImageUrl,
        width: 1200,
        height: 630,
        type: "image/png",
        alt: socialImageAlt
      }
    ]
  },
  twitter: {
    card: "summary_large_image",
    title: "DotMatch - Clear Read Assignment",
    description:
      "See which sequencing reads match known DNA targets and which remain ambiguous, unmatched, or invalid.",
    images: [
      {
        url: twitterImageUrl,
        width: 1200,
        height: 630,
        alt: socialImageAlt
      }
    ]
  }
};

export const viewport: Viewport = {
  width: "device-width",
  initialScale: 1
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
