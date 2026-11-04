import { describe, expect, it } from "vitest";
import { EMPTY_DRAFT, fromDraft, parseConstants, parseGroups, toDraft } from "./customCompounds";

describe("custom compound form", () => {
  it("reads groups with parentheses and commas", () => {
    expect(parseGroups("-CH3: 1, -CH2-: 1, -OH (alcohol): 1")).toEqual({ "-CH3": 1, "-CH2-": 1, "-OH (alcohol)": 1 });
    expect(parseGroups("CH3=2\nCH2COO = 1")).toEqual({ CH3: 2, CH2COO: 1 });
    expect(() => parseGroups("CH3: 1.5")).toThrow(/whole number/);
  });

  it("keeps units on constants and plain numbers as numbers", () => {
    expect(parseConstants("Tb: 78.4 degC; Pc: 61.4 bar, omega: 0.64")).toEqual({ Tb: "78.4 degC", Pc: "61.4 bar", omega: 0.64 });
    expect(() => parseConstants("Tcrit: 500")).toThrow(/unknown constant/);
  });

  it("builds a definition or says what is missing", () => {
    const draft = { ...EMPTY_DRAFT, name: "my-ethanol", formula: "C2H6O", joback: "-CH3: 1, -CH2-: 1, -OH (alcohol): 1" };
    const c = fromDraft(draft, ["water"]);
    expect(c).toEqual({ name: "my-ethanol", formula: "C2H6O", joback_groups: { "-CH3": 1, "-CH2-": 1, "-OH (alcohol)": 1 } });
    expect(fromDraft({ ...draft, name: "Water" }, ["water"])).toMatch(/already/);
    expect(fromDraft({ ...draft, formula: "" }, [])).toMatch(/formula/);
    expect(fromDraft({ ...draft, joback: "" }, [])).toMatch(/Joback groups/);
    expect(fromDraft({ ...draft, joback: "", constants: "Tb: 351, Tc: 514, Pc: 6.1e6" }, [])).toMatchObject({
      constants: { Tb: 351, Tc: 514, Pc: 6.1e6 },
    });
    expect(toDraft(c as never).joback).toBe("-CH3: 1, -CH2-: 1, -OH (alcohol): 1");
  });
});
