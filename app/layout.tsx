import type { Metadata } from "next";
import "./globals.css";

const siteUrl = process.env.NEXT_PUBLIC_SITE_URL ?? "https://dnncha.github.io/dotmatch";
const socialImageUrl = `${siteUrl}/dotmatch-og.png`;
const twitterImageUrl = `${siteUrl}/dotmatch-twitter.png`;
const socialImageAlt =
  "DotMatch social preview showing fixed-window barcode and guide assignment with auditable outcomes";

export const metadata: Metadata = {
  metadataBase: new URL(siteUrl),
  applicationName: "DotMatch",
  title: "DotMatch - Barcode Autopsy for Fixed-Window FASTQs",
  description:
    "DotMatch turns known short-DNA FASTQs into barcode splits, CRISPR guide counts, and QC reports with explicit ambiguity and autopsy artifacts.",
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
    "barcode demultiplexing",
    "barcode autopsy",
    "edit distance",
    "known-target assignment"
  ],
  openGraph: {
    title: "DotMatch - Barcode Autopsy for Fixed-Window FASTQs",
    description:
      "Barcode splits, CRISPR guide counts, and QC reports for known short-DNA targets, with ambiguous reads reported instead of guessed.",
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
    title: "DotMatch - Barcode Autopsy for Fixed-Window FASTQs",
    description:
      "Barcode splits, CRISPR guide counts, and QC reports for known short-DNA targets, with ambiguous reads reported instead of guessed.",
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
