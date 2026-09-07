# Tekstūros požymių analizė ir klasifikavimas

Šiame projekte atliekama tekstūros požymių analizė, jų paruošimas ir atranka, 
klasifikavimo bei klasterizavimo metodų taikymas ir rezultatų palyginimas 
su medicininėmis metrikomis.

## Projekto struktūra

Projekte naudojami trys pagrindiniai R skriptai:

### 1. `Požymių_vertinimas_ir_klasifikavimas.R`

Skriptas skirtas:
- tekstūros požymių nuskaitymui;
- duomenų paruošimui ir pirminiam apdorojimui;
- požymių atrankai;
- klasifikavimo modelių taikymui;
- klasifikavimo modelių rezultatų vertinimui.

### 2. `Klasterizavimas_ir_jo_vertinimas_dvi_klasės.R`

Skriptas skirtas:
- tekstūros požymių projekcijų analizei;
- klasterizavimo metodų taikymui;
- dviejų klasių vizualizavimui;
- klasių / klasterių vizualinio atsiskyrimo vertinimui.

### 3. `Medicininių_metrikų_vertinimas_ir_klasterių_interpretacija.R`

Skriptas skirtas:
- bendram tiriamųjų klasterizavimui pagal tekstūros požymius;
- gautų klasterių interpretacijai;
- klasterių tarpusavio palyginimui;
- klasterių rezultatų palyginimui su medicininėmis metrikomis.

## Analizės eiga

```text
Tekstūros požymiai
        │
        ▼
Duomenų nuskaitymas ir paruošimas
        │
        ▼
Požymių atranka
        │
        ├──────────────► Klasifikavimo modeliai
        │                       │
        │                       ▼
        │                Modelių vertinimas
        │
        ▼
     Klasterizavimas
        │
        ├──────────────► Projekcijų analizė
        │
        ├──────────────► Dviejų klasių vizualinis
        │                 atsiskyrimo vertinimas
        │
        ▼
Klasterių interpretacija
        │
        ▼
Palyginimas su medicininėmis metrikomis

