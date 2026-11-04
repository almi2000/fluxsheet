import { useEffect, useRef, useState } from "react";
import type { RecentEntry } from "./persistence";

/** The toolbar's Open button: open a file, or one of the recently opened or saved flowsheets. */
export function OpenMenu({
  recent,
  onOpenFile,
  onOpenRecent,
  onClearRecent,
}: {
  recent: RecentEntry[];
  onOpenFile: () => void;
  onOpenRecent: (entry: RecentEntry) => void;
  onClearRecent: () => void;
}) {
  const [open, setOpen] = useState(false);
  const ref = useRef<HTMLDivElement>(null);
  useEffect(() => {
    if (!open) return;
    const close = (e: MouseEvent) => {
      if (!ref.current?.contains(e.target as Node)) setOpen(false);
    };
    const esc = (e: KeyboardEvent) => e.key === "Escape" && setOpen(false);
    document.addEventListener("mousedown", close);
    document.addEventListener("keydown", esc);
    return () => {
      document.removeEventListener("mousedown", close);
      document.removeEventListener("keydown", esc);
    };
  }, [open]);

  const pick = (f: () => void) => () => {
    setOpen(false);
    f();
  };
  return (
    <div className="open-menu" ref={ref}>
      <button onClick={() => setOpen((o) => !o)} aria-haspopup="menu" aria-expanded={open}>
        Open… ▾
      </button>
      {open && (
        <div className="menu" role="menu">
          <button role="menuitem" onClick={pick(onOpenFile)}>Open file…</button>
          <div className="menu-label">Recent</div>
          {recent.length === 0 && <div className="menu-empty muted">Flowsheets you open or save appear here.</div>}
          {recent.map((e) => (
            <button
              key={`${e.name}-${e.at}`}
              role="menuitem"
              className="recent"
              title="A copy kept in this browser when the flowsheet was opened or saved"
              onClick={pick(() => onOpenRecent(e))}
            >
              <span>{e.name}</span>
              <small className="muted">{new Date(e.at).toLocaleString()}</small>
            </button>
          ))}
          {recent.length > 0 && (
            <button role="menuitem" className="menu-clear" onClick={pick(onClearRecent)}>
              Clear recent
            </button>
          )}
        </div>
      )}
    </div>
  );
}
