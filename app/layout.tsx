import type { Metadata, Viewport } from "next";
import "./globals.css";

const siteUrl = process.env.NEXT_PUBLIC_SITE_URL ?? "https://dnncha.github.io/dotmatch";
const socialImageUrl = `${siteUrl}/dotmatch-og.png`;
const twitterImageUrl = `${siteUrl}/dotmatch-twitter.png`;
const socialImageAlt =
  "DotMatch workflow showing known-target assignment outcomes from FASTQ reads";

export const metadata: Metadata = {
  metadataBase: new URL(siteUrl),
  applicationName: "DotMatch",
  title: "DotMatch - Known-target assignment from FASTQ",
  description:
    "Assign fixed FASTQ read windows to known guides, barcodes, feature tags, primers, and other short DNA targets without hiding ambiguous reads.",
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
    "known-target assignment",
    "barcode demultiplexing",
    "barcode panel design",
    "barcode troubleshooting",
    "FASTQ sequence matching",
    "CRISPR guide counting",
    "feature barcode assignment"
  ],
  openGraph: {
    title: "DotMatch - Known-target assignment from FASTQ",
    description:
      "Count guides, split inline barcodes, and inspect ambiguous or unmatched reads with an explicit target list.",
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
    title: "DotMatch - Known-target assignment from FASTQ",
    description:
      "Count guides, split inline barcodes, and inspect ambiguous or unmatched reads with an explicit target list.",
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
