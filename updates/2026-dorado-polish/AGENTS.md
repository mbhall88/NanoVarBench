# Dorado polish vs Clair3 (NanoVarBench update)

See `PLAN.md` for the agreed plan. Glossary: `CONTEXT.md`. Decisions: `docs/adr/`.

## Agent skills

### Issue tracker

GitHub Issues on `mbhall88/NanoVarBench` via the `gh` CLI (`gh` infers the repo inside this clone; `-R mbhall88/NanoVarBench` still works from elsewhere). See `docs/agents/issue-tracker.md`.

### Triage labels

Default five-role vocabulary (`needs-triage`, `needs-info`, `ready-for-agent`, `ready-for-human`, `wontfix`). See `docs/agents/triage-labels.md`.

### Domain docs

Single-context: `CONTEXT.md` + `docs/adr/` at this project's root. See `docs/agents/domain.md`.
