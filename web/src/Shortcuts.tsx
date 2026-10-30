import { useEffect, useRef } from "react";

const mod = typeof navigator !== "undefined" && /Mac|iPhone|iPad/.test(navigator.platform) ? "⌘" : "Ctrl";

export const SHORTCUTS: [keys: string[], what: string][] = [
  [[mod, "Z"], "Undo"],
  [[mod, "Shift", "Z"], "Redo (also Ctrl Y)"],
  [[mod, "C"], "Copy selected units"],
  [[mod, "X"], "Cut selected units"],
  [[mod, "V"], "Paste (at the mouse pointer)"],
  [[mod, "D"], "Duplicate selected units"],
  [[mod, "A"], "Select all units"],
  [["Del"], "Delete selection"],
  [["Shift", "click"], "Add a unit to the selection"],
  [["Shift", "drag"], "Box selection"],
  [["Esc"], "Clear selection"],
  [["?"], "Show this list"],
];

/** Popover listing the editor's keyboard shortcuts. */
export function ShortcutHelp({ onClose }: { onClose: () => void }) {
  const ref = useRef<HTMLDivElement>(null);
  useEffect(() => {
    const close = (e: MouseEvent) => {
      const t = e.target as HTMLElement;
      if (!ref.current?.contains(t) && !t.closest("[data-shortcut-toggle]")) onClose();
    };
    document.addEventListener("mousedown", close);
    return () => document.removeEventListener("mousedown", close);
  }, [onClose]);
  return (
    <div className="shortcut-help" ref={ref} role="dialog" aria-label="Keyboard shortcuts">
      <h2>Keyboard shortcuts</h2>
      <table>
        <tbody>
          {SHORTCUTS.map(([keys, what]) => (
            <tr key={what}>
              <td>
                {keys.map((k, i) => (
                  <span key={k}>
                    {i > 0 && " "}
                    <kbd>{k}</kbd>
                  </span>
                ))}
              </td>
              <td>{what}</td>
            </tr>
          ))}
        </tbody>
      </table>
      <p className="muted small">Shortcuts are off while you type in a field.</p>
    </div>
  );
}
