#ifndef COMREADER_HH
#define COMREADER_HH
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <set>
#include <bitset>
#include <unordered_map>
#include "TH1D.h"
#include "TTree.h"
#include "TFile.h"
#include "TString.h"

enum class WorkMode{
    WAVE,
    AMP,
    TEMP,
    MIX
};

enum class DetType{
    CALO,
    ITK
};

struct PairHash {
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);
        
        // 简单的哈希组合方法
        return h1 ^ h2;
        
        // 更好的哈希组合方法（来自Boost）
        // return h1 ^ (h2 << 1);
    }
};
class ComReader{
public:
static ComReader& getInstance() {
        // C++11 起，局部静态变量在首次调用时初始化，并自动保证线程安全
        static ComReader instance;
        return instance;
    }
    ComReader(const ComReader&) = delete;
    ComReader& operator=(const ComReader&) = delete;
    //Decode
    void decode(const std::string& filename);
    void decodeMix(const std::string& filename);
private:
    // 私有构造函数，防止外部实例化
    ComReader()=default;
    ~ComReader();

    //Logger
    std::shared_ptr<spdlog::logger> logger = spdlog::stdout_color_mt("ComReader");
    //IO
    bool openFile(const std::string& filename);
    std::ifstream *file;
    TFile *fout;
    std::string oname; //output file name
    void getOutputName(const std::string& prefix,const std::string& filename);

    //Initialize output tree for Calo, CsITK and HouseKeeping
    void initCaloTree();
    TTree *caloTree;

    void initCsITree();
    TTree *csiTree;

    void initHKTree(); //HouseKeeping Tree
    TTree *hkTree;

    //Find head
    bool findHead();
    bool readFEE();
    bool readCalo();
    bool readCaloWave();
    bool readBuffer(char *b,const int& NCHN);
    bool readCsI();
    bool readCsIWave();
	int readAmpTriggerIDMM();
	int readWaveTriggerIDMM();
    bool readHK();


    //Package attributes
    int PackageLength; //Package length
    int PackageID=-1;
    int TriggerID=-1;
    int TriggerID_csi=-1;
	// add Micromegas
	int TriggerIDMM = -1;
	int TriggerIDMM_csi = -1;

    int FEEID=-1;
    int nCaloHead=0;
    int nCsIHead=0;
    int nHKHead=0;
    int nCaliHead=0;
    int nhead=0;

    int Calo_EventID=0;
    int CsI_EventID=0;
    int Calo_EventCount=0;

    //Scientific Data
    //Calo
    std::vector<int> CellID; // CryID-FEEID-MBID-GID-CHNID
    std::vector<int> CellADC;
    std::vector<int> CellPLAT;
    //CsI

    //HouseKeeping Data
    int TPointID=0;
    int HK_EventCount=0;
    std::vector<float> C0;
    std::vector<float> C1;
    std::vector<float> C2;
    std::vector<float> T0;
    std::vector<float> T1;
    std::vector<float> T2;
    std::vector<float> T3;

    //Calo data buffer
    std::vector<std::vector<uint8_t>> calo_buffer;

    //Clear function
    void clear();
    void clearHK();

    bool findTempHead();
    int findMixHead(); // return 0 for scientific, 1 for temperature
    bool readWave();
    bool readAmp();
    bool readTemp();
    bool readMixTemp();
    bool readAmpCsI();
    bool readWaveCsI();

    //Get output name
    
    //Input file and output TFile and TTree


    //Data
    // static constexpr int NChn=26;
    // float chn[NChn];
    // float chnplat[NChn];

    // float npackages=0.;
    // int nhead=0;
    // int nevents_skipped=0;

    // float MixTemperature; // Temperature for MixReader
    // float MixCurrent; // Current for MixReader

    // int MixTempEventCount=0;
    // int MixTime;
    // std::set<int> set_FEEID;

    // //Temperature Current Data
    
    // static constexpr int NFEE = 4;

    // int MixFEEID=-1;
    

    

    //Constants
    static constexpr double e=1.602176634e-19;
    static constexpr double lowgain=0.549; //ADC/fC
    static constexpr double highgain=5.614; //ADC/fC
    static constexpr double MIP_Edep=36.; //MeV
    static constexpr double gain=50.;

    std::unordered_map<std::pair<int, int>, std::pair<int, int>, PairHash> channelMap_Calo = { //<FEEID,ChnID> <GID,CrystalID>
    {{3, 0}, {1, 1}}, {{3, 13}, {0, 1}}, {{4, 0}, {1, 1}}, {{4, 13}, {0, 1}},
    {{3, 1}, {1, 2}}, {{3, 14}, {0, 2}}, {{4, 1}, {1, 2}}, {{4, 14}, {0, 2}},
    {{1, 0}, {1, 3}},  {{1, 13}, {0, 3}}, {{2, 0}, {1, 3}},  {{2, 13}, {0, 3}},
    {{1, 1}, {1, 4}},  {{1, 14}, {0, 4}},  {{2, 1}, {1, 4}},  {{2, 14}, {0, 4}},
    {{1, 2}, {1, 5}},  {{1, 15}, {0, 5}},  {{2, 2}, {1, 5}},  {{2, 15}, {0, 5}},
    {{3, 2}, {1, 6}}, {{3, 15}, {0, 6}}, {{4, 2}, {1, 6}}, {{4, 15}, {0, 6}},
    {{3, 3}, {1, 7}}, {{3, 16}, {0, 7}}, {{4, 3}, {1, 7}}, {{4, 16}, {0, 7}},
    {{3, 4}, {1, 8}},  {{3, 17}, {0, 8}}, {{4, 4}, {1, 8}},  {{4, 17}, {0, 8}},
    {{1, 3}, {1, 9}},  {{1, 16}, {0, 9}}, {{2, 3}, {1, 9}},  {{2, 16}, {0, 9}},
    {{1, 4}, {1, 10}}, {{1, 17}, {0, 10}}, {{2, 4}, {1, 10}}, {{2, 17}, {0, 10}},
    {{3, 5}, {1, 11}},{{3, 18}, {0, 11}}, {{4, 5}, {1, 11}},{{4, 18}, {0, 11}},
    {{3, 6}, {1, 12}},{{3, 19}, {0, 12}},{{4, 6}, {1, 12}},{{4, 19}, {0, 12}},
    {{1, 5}, {1, 13}},{{1, 18}, {0, 13}},{{2, 5}, {1, 13}},{{2, 18}, {0, 13}},
    {{1, 6}, {1, 14}},{{1, 19}, {0, 14}},{{2, 6}, {1, 14}},{{2, 19}, {0, 14}},
    {{1, 7}, {1, 15}},{{1, 20}, {0, 15}}, {{2, 7}, {1, 15}},{{2, 20}, {0, 15}},
    {{3, 7}, {1, 16}},{{3, 20}, {0, 16}}, {{4, 7}, {1, 16}},{{4, 20}, {0, 16}},
    {{3, 8}, {1, 17}}, {{3, 21}, {0, 17}},{{4, 8}, {1, 17}}, {{4, 21}, {0, 17}},
    {{3, 9}, {1, 18}}, {{3, 22}, {0, 18}},{{4, 9}, {1, 18}}, {{4, 22}, {0, 18}},
    {{1, 8}, {1, 19}},{{1, 21}, {0, 19}},{{2, 8}, {1, 19}},{{2, 21}, {0, 19}},
    {{1, 9}, {1, 20}},{{1, 22}, {0, 20}}, {{2, 9}, {1, 20}},{{2, 22}, {0, 20}},
    {{3, 10}, {1, 21}}, {{3, 23}, {0, 21}}, {{4, 10}, {1, 21}}, {{4, 23}, {0, 21}},
    {{3, 11}, {1, 22}}, {{3, 24}, {0, 22}}, {{4, 11}, {1, 22}}, {{4, 24}, {0, 22}},
    {{1, 10}, {1, 23}},{{1, 23}, {0, 23}},{{2, 10}, {1, 23}},{{2, 23}, {0, 23}},
    {{1, 11}, {1, 24}},{{1, 24}, {0, 24}},{{2, 11}, {1, 24}},{{2, 24}, {0, 24}},
    {{1, 12}, {1, 25}},{{1, 25}, {0, 25}},{{2, 12}, {1, 25}},{{2, 25}, {0, 25}}
};

    std::unordered_map<std::pair<int, int>, std::pair<int, int>, PairHash> channelMap_CsI = { //<FEEID,ChnID> <GID,CrystalID>
    {{3, 22}, {1, 1}}, {{3, 25}, {0, 1}}, {{4, 22}, {1, 1}}, {{4, 25}, {0, 1}},
    {{3, 14}, {1, 2}}, {{3, 21}, {0, 2}}, {{4, 14}, {1, 2}}, {{4, 21}, {0, 2}},
    {{1, 4}, {1, 3}},  {{1, 17}, {0, 3}}, {{2, 4}, {1, 3}},  {{2, 17}, {0, 3}},
    {{1, 2}, {1, 4}},  {{1, 9}, {0, 4}},  {{2, 2}, {1, 4}},  {{2, 9}, {0, 4}},
    {{1, 0}, {1, 5}},  {{1, 1}, {0, 5}},  {{2, 0}, {1, 5}},  {{2, 1}, {0, 5}},
    {{3, 20}, {1, 6}}, {{3, 23}, {0, 6}}, {{4, 20}, {1, 6}}, {{4, 23}, {0, 6}},
    {{3, 12}, {1, 7}}, {{3, 19}, {0, 7}}, {{4, 12}, {1, 7}}, {{4, 19}, {0, 7}},
    {{3, 6}, {1, 8}},  {{3, 17}, {0, 8}}, {{4, 6}, {1, 8}},  {{4, 17}, {0, 8}},
    {{1, 8}, {1, 9}},  {{1, 11}, {0, 9}}, {{2, 8}, {1, 9}},  {{2, 11}, {0, 9}},
    {{1, 6}, {1, 10}}, {{1, 3}, {0, 10}}, {{2, 6}, {1, 10}}, {{2, 3}, {0, 10}},
    {{3, 18}, {1, 11}},{{3, 7}, {0, 11}}, {{4, 18}, {1, 11}},{{4, 7}, {0, 11}},
    {{3, 10}, {1, 12}},{{3, 13}, {0, 12}},{{4, 10}, {1, 12}},{{4, 13}, {0, 12}},
    {{1, 10}, {1, 13}},{{1, 19}, {0, 13}},{{2, 10}, {1, 13}},{{2, 19}, {0, 13}},
    {{1, 14}, {1, 14}},{{1, 13}, {0, 14}},{{2, 14}, {1, 14}},{{2, 13}, {0, 14}},
    {{1, 20}, {1, 15}},{{1, 5}, {0, 15}}, {{2, 20}, {1, 15}},{{2, 5}, {0, 15}},
    {{3, 16}, {1, 16}},{{3, 5}, {0, 16}}, {{4, 16}, {1, 16}},{{4, 5}, {0, 16}},
    {{3, 8}, {1, 17}}, {{3, 11}, {0, 17}},{{4, 8}, {1, 17}}, {{4, 11}, {0, 17}},
    {{3, 4}, {1, 18}}, {{3, 15}, {0, 18}},{{4, 4}, {1, 18}}, {{4, 15}, {0, 18}},
    {{1, 16}, {1, 19}},{{1, 15}, {0, 19}},{{2, 16}, {1, 19}},{{2, 15}, {0, 19}},
    {{1, 22}, {1, 20}},{{1, 7}, {0, 20}}, {{2, 22}, {1, 20}},{{2, 7}, {0, 20}},
    {{3, 2}, {1, 21}}, {{3, 3}, {0, 21}}, {{4, 2}, {1, 21}}, {{4, 3}, {0, 21}},
    {{3, 0}, {1, 22}}, {{3, 9}, {0, 22}}, {{4, 0}, {1, 22}}, {{4, 9}, {0, 22}},
    {{1, 12}, {1, 23}},{{1, 25}, {0, 23}},{{2, 12}, {1, 23}},{{2, 25}, {0, 23}},
    {{1, 18}, {1, 24}},{{1, 23}, {0, 24}},{{2, 18}, {1, 24}},{{2, 23}, {0, 24}},
    {{1, 24}, {1, 25}},{{1, 21}, {0, 25}},{{2, 24}, {1, 25}},{{2, 21}, {0, 25}}
};

};

#endif
