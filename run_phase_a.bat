@echo off
mkdir results 2>nul
for %%G in (INS ACTB MYH9 ALB GAPDH TNF VEGFA APOE CFTR MYC PTEN HSP90AA1) do (
    echo Running %%G...
    python pipeline.py %%G > results\%%G.txt 2>&1
)
echo Done! All results saved in the results\ folder.