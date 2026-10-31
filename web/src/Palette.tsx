import { UnitIcon } from "./icons";
import type { UnitType } from "./types";

export function Palette({ types }: { types: UnitType[] }) {
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
                  title={t.description}
                  onDragStart={(e) => {
                    e.dataTransfer.setData("application/fluxsheet-unit", t.type);
                    e.dataTransfer.effectAllowed = "move";
                  }}
                >
                  <UnitIcon type={t.type} size={38} />
                  <span>{t.label}</span>
                </div>
              ))}
          </div>
        </section>
      ))}
    </aside>
  );
}
