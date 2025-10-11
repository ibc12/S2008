#include "ActInputParser.h"
#include "ActModularDetector.h"
#include "ActSilDetector.h"

void checkActions(){
    ActRoot::InputParser in {"./configs/detector.conf"};
    ActRoot::ModularDetector mod;
    mod.ReadConfiguration(in.GetBlock("Modular"));
    // mod.GetParameters()->Print();
    ActRoot::SilDetector sil;
    sil.ReadConfiguration(in.GetBlock("Silicons"));
    sil.GetParameters()->Print();
}
