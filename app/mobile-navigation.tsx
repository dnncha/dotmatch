"use client";

import { useRef } from "react";

type MobileNavigationProps = {
  links: readonly (readonly [label: string, href: string])[];
  docsUrl: string;
  repoUrl: string;
};

export function MobileNavigation({ links, docsUrl, repoUrl }: MobileNavigationProps) {
  const detailsRef = useRef<HTMLDetailsElement>(null);

  function closeMenu() {
    detailsRef.current?.removeAttribute("open");
  }

  return (
    <details ref={detailsRef} className="mobile-nav">
      <summary>Menu</summary>
      <nav aria-label="Mobile navigation">
        {links.map(([label, href]) => (
          <a key={href} href={href} onClick={closeMenu}>{label}</a>
        ))}
        <a href={docsUrl} onClick={closeMenu}>Docs</a>
        <a href={repoUrl} onClick={closeMenu}>GitHub</a>
      </nav>
    </details>
  );
}
