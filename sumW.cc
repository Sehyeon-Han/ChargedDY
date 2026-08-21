#include <iostream>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

void sumW(const char* inFile,
          const char* outFile = "sumW.txt")
{
    TH1::SetDefaultSumw2();
    
    TChain chain("Events");
    chain.Add(inFile);

    TTreeReader reader(&chain);
    TTreeReaderValue<Float_t> genWeight(reader, "genWeight");

    double sumW = 0.0;

    while(reader.Next())
    {
        sumW += (double)(*genWeight);
    }

    FILE* f = fopen(outFile, "w");
    fprintf(f, "%f\n", sumW);
    fclose(f);
}