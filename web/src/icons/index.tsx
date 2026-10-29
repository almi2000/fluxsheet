// Unit icons. Every file in this folder that exports `shapes` (unit type -> SVG content on
// a 40x40 grid) is picked up automatically, so each unit family keeps its own icon file.
import { S } from "./style";

const modules = import.meta.glob<{ shapes?: Record<string, JSX.Element> }>("./*.tsx", { eager: true });
const shapes: Record<string, JSX.Element> = Object.assign(
  {},
  ...Object.values(modules).map((m) => m.shapes ?? {}),
);

export function UnitIcon({ type, size = 40 }: { type: string; size?: number }) {
  return (
    <svg width={size} height={size} viewBox="0 0 40 40" aria-hidden="true">
      {shapes[type] ?? <rect {...S} x="6" y="6" width="28" height="28" rx="4" />}
    </svg>
  );
}
