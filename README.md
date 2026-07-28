# ChatGEM

Code and data files accompanying the chatGEM manuscript. 

chatGEM accelerates genome scale metabolic modeling workflows via a grounded coding agent 


## Repository layout

```
RAG_COBRApy_codes/       # 52 numbered COBRApy scripts + model files and data assets
ecGEM_ChatGEM/           # per-experiment enzyme-constrained GEM runs for P. putida KT2440 OPAL strains
iBAQ/                    # iBAQ proteomics processing scripts and raw data
/skills/chatgem/         # ChatGEM agent skill --> move to .claude to use with Claude code
```

## Setup

```bash
git clone https://github.com/pnnl/ChatGEM
cd ChatGEM
uv sync
```


## Related Work

ADEPT agentic platform and core system: [https://github.com/pnnl/adept-agentic](https://github.com/pnnl/adept-agentic) and [https://github.com/pnnl/adept-agentic-framework-core](https://github.com/pnnl/adept-agentic-framework-core)


ChatGEM preprint: [https://www.biorxiv.org/content/10.64898/2026.07.20.739662v2](https://www.biorxiv.org/content/10.64898/2026.07.20.739662v2)



---
Niaz Chowdhury & August George, 2026, PNNL