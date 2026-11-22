import { UnitIcon } from "./icons";
import type { UnitType } from "./types";

export function Palette({ types, problems = [] }: { types: UnitType[]; problems?: string[] }) {
  const categories = [...new Set(types.map((t) => t.category))];
  return (
    <aside className="palette">
      <h2>Units</h2>
      <p className="muted small">Drag onto the canvas</p>
      {categories.map((cat) => (
        <section key={cat}>
          <h3>{cat}</h3>
          <div className="palette-grid">
            {types
              .filter((t) => t.category === cat)
              .map((t) => (
                <div
                  key={t.type}
                  className="palette-item"
                  draggable
                  title={t.plugin ? `${t.description}

From the plug-in ${t.plugin}.` : t.description}
                  onDragStart={(e) => {
                    e.dataTransfer.setData("application/fluxsheet-unit", t.type);
                    e.dataTransfer.effectAllowed = "move";
                  }}
                >
                  <UnitIcon type={t.type} size={38} />
                  <span>{t.label}</span>
                  {t.plugin && <span className="palette-plugin">plug-in</span>}
                </div>
              ))}
          </div>
        </section>
      ))}
      {problems.length > 0 && (
        <section className="palette-problems">
          <h3>Plug-ins not loaded</h3>
          <ul>
            {problems.map((p) => (
              <li key={p}>{p}</li>
            ))}
          </ul>
          <p className="muted small">
            Run <code>fluxsheet plugins</code> for details.
          </p>
        </section>
      )}
    </aside>
  );
}
