import { renderToStaticMarkup } from "react-dom/server";
import { expect, it } from "vitest";
import { StreamInspector } from "./Inspector";
import type { StreamSettings } from "./types";

it("opens a recycle stream with a null start guess in a flowsheet file", () => {
  const settings: StreamSettings = JSON.parse('{"circulation":{"mass_flow":"0.1 kg/s"},"recycle_guess":null}');
  const html = renderToStaticMarkup(<StreamInspector
    components={["water"]} editable id="RETURN" onChange={() => {}} settings={settings}
  />);
  expect(html).toContain("Circulating mass flow");
  expect(html).toContain("Recycle starting values");
  expect(html).toContain('value="0.1"');
});
