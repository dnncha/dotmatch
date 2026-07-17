import type { Metadata, Viewport } from "next";
import "./globals.css";

const siteUrl = process.env.NEXT_PUBLIC_SITE_URL ?? "https://dnncha.github.io/dotmatch";
const socialImageUrl = `${siteUrl}/dotmatch-og.png`;
const twitterImageUrl = `${siteUrl}/dotmatch-twitter.png`;
const socialImageAlt =
  "AssayCode preview showing the DotMatch assignment outcomes for known-target sequencing assays";

export const metadata: Metadata = {
  metadataBase: new URL(siteUrl),
  applicationName: "AssayCode",
  title: "AssayCode - Design, Decode, and Diagnose Known-Target Assays",
  description:
    "AssayCode compiles, validates, decodes, and diagnoses known-target sequencing assays with the DotMatch engine.",
  authors: [{ name: "DotMatch maintainers", url: "https://github.com/dnncha/dotmatch" }],
  creator: "DotMatch maintainers",
  publisher: "AssayCode",
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
    "assay compiler",
    "AssayScript",
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
    title: "AssayCode - Assay Reliability",
    description:
      "Design the assay and trust the assignment with AssayCode, powered by DotMatch.",
    type: "website",
    siteName: "AssayCode",
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
    title: "AssayCode - Assay Reliability",
    description:
      "Design the assay and trust the assignment with AssayCode, powered by DotMatch.",
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
