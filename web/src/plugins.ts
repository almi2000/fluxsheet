import type { Method, PluginInfo } from "./types";

/** Option text of a property method; a plug-in's method says so. */
export function methodLabel(m: Method): string {
  return m.plugin ? `${m.label} (plug-in)` : m.label;
}

/** One line per plug-in that failed to load, for the palette. */
export function pluginProblems(plugins: PluginInfo[]): string[] {
  return plugins
    .filter((p) => p.error)
    .map((p) => `${p.distribution ?? p.value} (${p.name}): ${p.error}`);
}
