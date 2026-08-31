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
  title: "DotMatch — Known-target sequencing read assignment",
  description:
    "Match a fixed window in each sequencing read to known short DNA targets while keeping ambiguous, unmatched, and invalid reads visible.",
  authors: [{ name: "DotMatch maintainers", url: "https://github.com/dnncha/dotmatch" }],
  creator: "DotMatch maintainers",
  publisher: "DotMatch",
  category: "Bioinformatics software",
  alternates: {
    canonical: siteUrl
  },
  keywords: [
    "bioinformatics",
    "CRISPR",
    "FASTQ",
    "known-target sequencing",
    "barcode demultiplexing",
    "feature-barcode assignment",
    "guide capture",
    "Perturb-seq",
    "barcode panel design",
    "FASTQ sequence matching"
  ],
  openGraph: {
    title: "DotMatch — Match reads without hiding uncertainty",
    description:
      "Assign fixed read windows to known DNA targets while keeping ambiguous, unmatched, and invalid reads visible.",
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
    title: "DotMatch — Match reads without hiding uncertainty",
    description:
      "Assign fixed read windows to known DNA targets while keeping ambiguous, unmatched, and invalid reads visible.",
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
      <head>
        <link rel="describedby" href={`${siteUrl}/llms.txt`} />
      </head>
      <body>{children}</body>
    </html>
  );
}
