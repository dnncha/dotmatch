import type { Metadata, Viewport } from "next";
import "./globals.css";

const siteUrl = process.env.NEXT_PUBLIC_SITE_URL ?? "https://dnncha.github.io/dotmatch";
const socialImageUrl = `${siteUrl}/dotmatch-og.png`;
const twitterImageUrl = `${siteUrl}/dotmatch-twitter.png`;
const socialImageAlt =
  "DotMatch social preview showing CRISPR guide counts, barcode panel checks, barcode splits, and visible QC outcomes";

export const metadata: Metadata = {
  metadataBase: new URL(siteUrl),
  applicationName: "DotMatch",
  title: "DotMatch - Barcode Panels, Guide Counts, and Barcode QC",
  description:
    "DotMatch designs barcode panels, counts CRISPR guides, and splits fixed-position barcodes from FASTQ reads, with ambiguity kept visible.",
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
    "barcode panel design",
    "barcode troubleshooting",
    "edit distance",
    "known-target assignment"
  ],
  openGraph: {
    title: "DotMatch - Barcode Panels, Guide Counts, and Barcode QC",
    description:
      "Barcode panel certificates, CRISPR guide counts, barcode splits, and QC reports for known short-DNA targets.",
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
    title: "DotMatch - Barcode Panels, Guide Counts, and Barcode QC",
    description:
      "Barcode panel certificates, CRISPR guide counts, barcode splits, and QC reports for known short-DNA targets.",
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
