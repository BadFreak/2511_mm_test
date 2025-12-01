#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <array>
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "TKey.h"
#include "TH2D.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TLegend.h"
#include "TMath.h"
#include "TString.h"
#include "TLegend.h"
#include "TF1.h"
#include "TText.h"

// Landau ⊗ Gaussian helper (aka Langau).
Double_t LandauGauss(Double_t *x, Double_t *par) {
    constexpr Double_t invSqrt2Pi = 0.3989422804014327;
    const Double_t mpv = par[0];
    const Double_t landauWidth = par[1];
    const Double_t gaussSigma = par[2];
    const Double_t amplitude = par[3];
    if (landauWidth <= 0 || gaussSigma <= 0) {
        return 0.;
    }
    const Double_t xVal = x[0];
    const Double_t rangeLow = xVal - 5.0 * gaussSigma;
    const Double_t rangeHigh = xVal + 5.0 * gaussSigma;
    constexpr int nSteps = 100;
    const Double_t step = (rangeHigh - rangeLow) / nSteps;
    Double_t sum = 0.;
    for (int i = 0; i < nSteps; ++i) {
        const Double_t xx = rangeLow + (i + 0.5) * step;
        const Double_t landau = TMath::Landau(xx, mpv, landauWidth, true);
        const Double_t arg = (xVal - xx) / gaussSigma;
        const Double_t gauss = TMath::Exp(-0.5 * arg * arg);
        sum += landau * gauss;
    }
    return amplitude * step * invSqrt2Pi / gaussSigma * sum;
}

void ana(std::string option) {
    struct CsiData {
        std::vector<int> cellADC;
        std::vector<int> cellPLAT;
    };
    struct AdasData {
        Double_t xSlope;
        Double_t xIntercept;
        Double_t ySlope;
        Double_t yIntercept;
    };
    
    const std::string outputDir = "result/result" + option;
    if (gSystem->AccessPathName(outputDir.c_str())) {
        if (gSystem->mkdir(outputDir.c_str(), true) != 0) {
            std::cerr << "错误：无法创建输出目录 " << outputDir << std::endl;
            return;
        }
    }
    
    // 数组1: 存储 csiTree 的 TriggerIDMM
    std::vector<int> triggerIDMMArray;
    std::map<int, std::vector<CsiData>> csiDataByTrigger;
    
    // 读取 result_decode.root 文件中的 csiTree 的 TriggerIDMM
    std::string datadecodeFile = outputDir + "/result_decode" + option + ".root";
    
    std::cout << "读取文件: " << datadecodeFile << std::endl;
    
    TFile *file = TFile::Open(datadecodeFile.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "错误：无法打开文件 " << datadecodeFile << std::endl;
        return;
    }
    
    // 直接获取 csiTree
    TTree *csiTree = (TTree*)file->Get("csiTree");
    if (!csiTree) {
        std::cerr << "错误：无法找到 csiTree" << std::endl;
        file->Close();
        delete file;
        return;
    }
    
    // 检查需要的 branch
    TBranch *branchTrigger = csiTree->GetBranch("TriggerIDMM");
    TBranch *branchCellADC = csiTree->GetBranch("CellADC");
    TBranch *branchCellPLAT = csiTree->GetBranch("CellPLAT");
    if (!branchTrigger || !branchCellADC || !branchCellPLAT) {
        std::cerr << "错误：csiTree 中缺少 TriggerIDMM/CellADC/CellPLAT branch" << std::endl;
        file->Close();
        delete file;
        return;
    }
    
    Int_t triggerIDMM;
    std::vector<int> *cellADCVec = nullptr;
    std::vector<int> *cellPLATVec = nullptr;
    csiTree->SetBranchAddress("TriggerIDMM", &triggerIDMM);
    csiTree->SetBranchAddress("CellADC", &cellADCVec);
    csiTree->SetBranchAddress("CellPLAT", &cellPLATVec);
    
    Long64_t nEntries = csiTree->GetEntries();
    std::cout << "  Tree: csiTree, 条目数: " << nEntries << std::endl;
    
    for (Long64_t i = 0; i < nEntries; i++) {
        csiTree->GetEntry(i);
        triggerIDMMArray.push_back(triggerIDMM);
        if (cellADCVec && cellPLATVec) {
            CsiData entry;
            entry.cellADC = *cellADCVec;
            entry.cellPLAT = *cellPLATVec;
            csiDataByTrigger[triggerIDMM].push_back(std::move(entry));
        }
		// cout triggerIDMM
		// std::cout << "TriggerID MM : " << triggerIDMM << std::endl;
    }
    
    file->Close();
    delete file;
    
    std::cout << "收集到 " << triggerIDMMArray.size() << " 个 TriggerIDMM" << std::endl;
    
    // 数组2: 存储 adas_track_data 文件中的 EventID
    std::vector<int> eventIDArray;
    std::map<int, std::vector<AdasData>> adasDataByEvent;
    
    // 直接读取指定的文件和 Tree
    std::string adasFile = outputDir + "/adas_track_data" + option + ".root";
    
    std::cout << "\n读取文件: " << adasFile << std::endl;
    
    TFile *adasFilePtr = TFile::Open(adasFile.c_str(), "READ");
    if (!adasFilePtr || adasFilePtr->IsZombie()) {
        std::cerr << "错误：无法打开文件 " << adasFile << std::endl;
        return;
    }
    
    // 直接获取 adas_track_data Tree
    TTree *adasTree = (TTree*)adasFilePtr->Get("adas_track_data");
    if (!adasTree) {
        std::cerr << "错误：无法找到 adas_track_data Tree" << std::endl;
        adasFilePtr->Close();
        delete adasFilePtr;
        return;
    }
    
    auto findBranch = [&](const std::vector<std::string> &candidates, std::string &chosen) -> TBranch* {
        for (const auto &name : candidates) {
            TBranch *br = adasTree->GetBranch(name.c_str());
            if (br) {
                chosen = name;
                return br;
            }
        }
        return nullptr;
    };
    
    std::string eventBranchName, xSlopeBranchName, xInterceptBranchName, ySlopeBranchName, yInterceptBranchName;
    TBranch *branchEvent = findBranch({"event_id", "EventID", "TriggerID"}, eventBranchName);
    TBranch *branchXSlope = findBranch({"x_slope", "xSlope", "xslope"}, xSlopeBranchName);
    TBranch *branchXIntercept = findBranch({"x_intercept", "x_Intercept", "xIntercept"}, xInterceptBranchName);
    TBranch *branchYSlope = findBranch({"y_slope", "ySlope", "yslope"}, ySlopeBranchName);
    TBranch *branchYIntercept = findBranch({"y_intercept", "y_Intercept", "yIntercept"}, yInterceptBranchName);
    
    if (!branchEvent || !branchXSlope || !branchXIntercept || !branchYSlope || !branchYIntercept) {
        std::cerr << "错误：adas_track_data Tree 中缺少必要 branch（event_id/x_slope/x_intercept/y_slope/y_intercept 的某个变体）" << std::endl;
        adasFilePtr->Close();
        delete adasFilePtr;
        return;
    }
    
    Int_t eventID;
    Float_t xSlopeF = 0, xInterceptF = 0, ySlopeF = 0, yInterceptF = 0;
    adasTree->SetBranchAddress(eventBranchName.c_str(), &eventID);
    adasTree->SetBranchAddress(xSlopeBranchName.c_str(), &xSlopeF);
    adasTree->SetBranchAddress(xInterceptBranchName.c_str(), &xInterceptF);
    adasTree->SetBranchAddress(ySlopeBranchName.c_str(), &ySlopeF);
    adasTree->SetBranchAddress(yInterceptBranchName.c_str(), &yInterceptF);
    
    Long64_t nEntriesEvent = adasTree->GetEntries();
    std::cout << "  Tree: adas_track_data, Branch: event_id, 条目数: " << nEntriesEvent << std::endl;
    
    for (Long64_t i = 0; i < nEntriesEvent; i++) {
        adasTree->GetEntry(i);
        eventIDArray.push_back(eventID);
        adasDataByEvent[eventID].push_back(
            {static_cast<Double_t>(xSlopeF), static_cast<Double_t>(xInterceptF),
             static_cast<Double_t>(ySlopeF), static_cast<Double_t>(yInterceptF)});
		// cout eventID
		// std::cout << "eventID : " << eventID << std::endl;
    }
    
    adasFilePtr->Close();
    delete adasFilePtr;
    
    std::cout << "收集到 " << eventIDArray.size() << " 个 EventID" << std::endl;
    
    // 找出公共数字
    // 将数组转换为 set 以便快速查找
    std::set<int> triggerIDMMSet(triggerIDMMArray.begin(), triggerIDMMArray.end());
    std::set<int> eventIDSet(eventIDArray.begin(), eventIDArray.end());
    
    // 找出公共数字
    std::vector<int> commonNumbers;
    for (int id : triggerIDMMSet) {
        if (eventIDSet.find(id) != eventIDSet.end()) {
            commonNumbers.push_back(id);
        }
    }
    
    // 排序公共数字
    std::sort(commonNumbers.begin(), commonNumbers.end());
    
    // 输出结果
    std::cout << "\n==========================================" << std::endl;
    std::cout << "统计结果：" << std::endl;
    std::cout << "==========================================" << std::endl;
    std::cout << "TriggerIDMM 数组大小: " << triggerIDMMArray.size() << std::endl;
    std::cout << "EventID 数组大小: " << eventIDArray.size() << std::endl;
    std::cout << "公共数字数量: " << commonNumbers.size() << std::endl;
    std::cout << "==========================================" << std::endl;
    
    // 输出文件保存公共数据
    if (!commonNumbers.empty()) {
        std::string outputFileName = outputDir + "/common_trigger_data_" + option + ".root";
        TFile *outputFile = TFile::Open(outputFileName.c_str(), "RECREATE");
        if (!outputFile || outputFile->IsZombie()) {
            std::cerr << "错误：无法创建输出文件 " << outputFileName << std::endl;
            delete outputFile;
            return;
        }
        
        TTree *commonTree = new TTree("commonTree", "Matched trigger data");
        Int_t outTriggerID = 0;
        Double_t outXSlope = 0, outXIntercept = 0, outYSlope = 0, outYIntercept = 0;
        Double_t outXPos = 0, outYPos = 0;
        std::vector<int> outCellADC;
        std::vector<int> outCellPLAT;
        std::vector<int> outCellRealADC;
        constexpr int gridDivision = 10;
        constexpr int cornerBlock = gridDivision / 2;
        // constexpr double gridXMin = 155.0;
        // constexpr double gridXMax = 455.0;
        constexpr double gridXMin = 160.0;
		constexpr double gridXMax = 460.0;
		constexpr double gridYMin = 140.5;
        constexpr double gridYMax = 440.5;
        constexpr double gridXSpan = gridXMax - gridXMin;
        constexpr double gridYSpan = gridYMax - gridYMin;
        constexpr double gridCellSizeX = gridXSpan / gridDivision;
        constexpr double gridCellSizeY = gridYSpan / gridDivision;
        TH2D *xyDistribution =
            new TH2D("XYDistribution", "Reconstructed XY;X position (mm);Y position (mm)",
                     gridDivision * 5, gridXMin, gridXMax, gridDivision * 5, gridYMin, gridYMax);
        xyDistribution->SetStats(0);
        std::vector<std::vector<TH1D*>> gridCellHistograms(
            gridDivision * gridDivision, std::vector<TH1D*>(8, nullptr));
        for (int padIndex = 0; padIndex < gridDivision * gridDivision; ++padIndex) {
            for (int ch = 0; ch < 8; ++ch) {
                std::string histName =
                    "CellRealADC_pad" + std::to_string(padIndex) + "_ch" + std::to_string(ch);
                std::string histTitle =
                    "Pad " + std::to_string(padIndex) + " channel " + std::to_string(ch);
                gridCellHistograms[padIndex][ch] =
                    new TH1D(histName.c_str(), histTitle.c_str(), 200, 200., 3000.);
                gridCellHistograms[padIndex][ch]->Sumw2();
            }
        }
        
        commonTree->Branch("TriggerIDMM", &outTriggerID, "TriggerIDMM/I");
        commonTree->Branch("x_slope", &outXSlope, "x_slope/D");
        commonTree->Branch("x_intercept", &outXIntercept, "x_intercept/D");
        commonTree->Branch("y_slope", &outYSlope, "y_slope/D");
        commonTree->Branch("y_intercept", &outYIntercept, "y_intercept/D");
        commonTree->Branch("x_pos", &outXPos, "x_pos/D");
        commonTree->Branch("y_pos", &outYPos, "y_pos/D");
        commonTree->Branch("CellADC", &outCellADC);
        commonTree->Branch("CellPLAT", &outCellPLAT);
        commonTree->Branch("CellRealADC", &outCellRealADC);
        
        size_t savedEntries = 0;
        for (int id : commonNumbers) {
            const auto &adasVec = adasDataByEvent[id];
            const auto &csiVec = csiDataByTrigger[id];
            for (const auto &adasEntry : adasVec) {
                for (const auto &csiEntry : csiVec) {
                    outTriggerID = id;
                    outXSlope = adasEntry.xSlope;
                    outXIntercept = adasEntry.xIntercept;
                    outYSlope = adasEntry.ySlope;
                    outYIntercept = adasEntry.yIntercept;
                    outXPos = outXSlope * 674.4 + outXIntercept;
                    outYPos = outYSlope * 674.4 + outYIntercept;
                    if (xyDistribution) {
                        xyDistribution->Fill(outXPos, outYPos);
                    }
                    outCellADC = csiEntry.cellADC;
                    outCellPLAT = csiEntry.cellPLAT;
                    outCellRealADC.clear();
                    const size_t nCells = std::min(outCellADC.size(), outCellPLAT.size());
                    outCellRealADC.reserve(nCells);
                    bool withinGrid = (outXPos >= gridXMin && outXPos < gridXMax &&
                                       outYPos >= gridYMin && outYPos < gridYMax);
                    int padIndex = -1;
                    if (withinGrid) {
                        int colFromRight =
                            static_cast<int>((outXPos - gridXMin) / gridCellSizeX);
                        int rowFromTop =
                            static_cast<int>((outYPos - gridYMin) / gridCellSizeY);
                        colFromRight = std::min(std::max(colFromRight, 0), gridDivision - 1);
                        rowFromTop = std::min(std::max(rowFromTop, 0), gridDivision - 1);
                        int colForCanvas = (gridDivision - 1) - colFromRight;
                        padIndex = rowFromTop * gridDivision + colForCanvas;
                    }
                    for (size_t idx = 0; idx < nCells; ++idx) {
                        int realADC = outCellADC[idx] - outCellPLAT[idx];
                        outCellRealADC.push_back(realADC);
                        if (withinGrid && padIndex >= 0 && idx < gridCellHistograms[padIndex].size()) {
                            if (realADC>300) gridCellHistograms[padIndex][idx]->Fill(realADC);
                        }
                    }
                    commonTree->Fill();
                    ++savedEntries;
                }
            }
        }
        
        commonTree->Write();
        
        std::array<int, 8> colorPalette = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1,
                                           kCyan + 1, kOrange + 7, kGray + 1, kBlack};
        
        // 对每个 gridCellHistograms 进行 Landau-Gaussian 拟合
        std::vector<std::vector<TF1*>> fitFunctions(
            gridDivision * gridDivision, std::vector<TF1*>(8, nullptr));
        // 存储 MPV 值：padIndex -> channel -> MPV
        std::vector<std::vector<Double_t>> mpvValues(
            gridDivision * gridDivision, std::vector<Double_t>(8, 0.0));
        for (int padIndex = 0; padIndex < gridDivision * gridDivision; ++padIndex) {
            for (int ch = 0; ch < 8; ++ch) {
                TH1D *hist = gridCellHistograms[padIndex][ch];
                if (!hist || hist->GetEntries() < 300) {
                    continue;  // 只对 entries >= 300 的进行拟合
                }
                
                // 获取 histogram 的最大值位置作为初始 MPV
                Int_t maxBin = hist->GetMaximumBin();
                Double_t maxX = hist->GetXaxis()->GetBinCenter(maxBin);
                Double_t maxY = hist->GetBinContent(maxBin);
                
                // 估算初始参数
                Double_t mpv = maxX;
                Double_t landauWidth = 100.0;  // 初始 Landau 宽度
                Double_t gaussSigma = 50.0;     // 初始高斯 sigma
                Double_t amplitude = maxY * 50.0;  // 初始幅度
                
                // 创建拟合函数
                std::string fitName = "fit_pad" + std::to_string(padIndex) + "_ch" + std::to_string(ch);
                TF1 *fitFunc = new TF1(fitName.c_str(), LandauGauss, 200., 3000., 4);
                fitFunc->SetParameters(mpv, landauWidth, gaussSigma, amplitude);
                fitFunc->SetParNames("MPV", "LandauWidth", "GaussSigma", "Amplitude");
                
                // 设置参数范围
                fitFunc->SetParLimits(0, 200., 3000.);      // MPV
                fitFunc->SetParLimits(1, 10., 500.);       // LandauWidth
                fitFunc->SetParLimits(2, 10., 200.);       // GaussSigma
                fitFunc->SetParLimits(3, 0., amplitude * 10.);  // Amplitude
                
                // 使用 histogram 本身的颜色
                Int_t histColor = colorPalette[ch % colorPalette.size()];
                fitFunc->SetLineColor(histColor);
                fitFunc->SetLineWidth(2);
                fitFunc->SetLineStyle(2);  // 虚线
                
                // 执行拟合
                Int_t fitStatus = hist->Fit(fitFunc, "QN0", "", 200., 3000.);
                if (fitStatus == 0) {
                    fitFunctions[padIndex][ch] = fitFunc;
                    // 获取拟合得到的 MPV 值
                    mpvValues[padIndex][ch] = fitFunc->GetParameter(0);
                } else {
                    delete fitFunc;
                }
            }
        }
        
        // 创建 MPV 统计的 canvas 和 TH2D
        struct MPVCornerInfo {
            std::string name;
            TCanvas *canvas = nullptr;
            std::vector<TH2D*> th2dChannels;  // 2 个 TH2D，每个对应一个 channel
            std::vector<int> channels;  // 该 corner 对应的 channel 编号
            int rowStart;
            int colStart;
        };
        std::vector<MPVCornerInfo> mpvCornerCanvases;
        
        // 定义 4 个 corner 的位置和对应的 channel
        struct CornerDef {
            std::string name;
            int rowStart;
            int colStart;
            std::vector<int> channels;
        };
        std::vector<CornerDef> cornerDefs = {
            {"MPV_TopLeft", 0, 0, {2, 3}},
            {"MPV_TopRight", 0, gridDivision - cornerBlock, {0, 1}},
            {"MPV_BottomLeft", gridDivision - cornerBlock, 0, {6, 7}},
            {"MPV_BottomRight", gridDivision - cornerBlock, gridDivision - cornerBlock, {4, 5}}
        };
        
        // 创建 4 个 corner canvas，每个有 2 个 pad（对应 2 个 channel）
        for (const auto &cornerDef : cornerDefs) {
            MPVCornerInfo cornerInfo;
            cornerInfo.name = cornerDef.name;
            cornerInfo.rowStart = cornerDef.rowStart;
            cornerInfo.colStart = cornerDef.colStart;
            cornerInfo.channels = cornerDef.channels;
            cornerInfo.canvas = new TCanvas(cornerInfo.name.c_str(), 
                (cornerInfo.name + " MPV Distribution").c_str(), 1600, 800);
            cornerInfo.canvas->Divide(2, 1);  // 2 列 1 行，横向排列，共 2 个 pad
            
            // 计算该 corner 的物理位置范围
            Double_t cornerXMin = gridXMin + cornerInfo.colStart * gridCellSizeX;
            Double_t cornerXMax = gridXMin + (cornerInfo.colStart + cornerBlock) * gridCellSizeX;
            Double_t cornerYMin = gridYMin + cornerInfo.rowStart * gridCellSizeY;
            Double_t cornerYMax = gridYMin + (cornerInfo.rowStart + cornerBlock) * gridCellSizeY;
            
            // 为每个 channel 创建 TH2D（使用物理位置坐标）
            for (int ch : cornerInfo.channels) {
                std::string th2dName = cornerInfo.name + "_ch" + std::to_string(ch);
                std::string th2dTitle = cornerInfo.name + " Channel " + std::to_string(ch) + " MPV;X position (mm);Y position (mm)";
                TH2D *th2d = new TH2D(th2dName.c_str(), th2dTitle.c_str(),
                    cornerBlock, cornerXMin, cornerXMax, cornerBlock, cornerYMin, cornerYMax);
                th2d->SetStats(0);
                cornerInfo.th2dChannels.push_back(th2d);
            }
            
            // 填充 TH2D：遍历该 corner 的所有 pad
            for (int dr = 0; dr < cornerBlock; ++dr) {
                for (int dc = 0; dc < cornerBlock; ++dc) {
                    int row = cornerInfo.rowStart + dr;
                    int col = cornerInfo.colStart + dc;
                    if (row < gridDivision && col < gridDivision) {
                        int padIdx = row * gridDivision + col;
                        // 计算该 pad 的物理位置（中心点）
                        Double_t padX = gridXMin + (col + 0.5) * gridCellSizeX;
                        Double_t padY = gridYMin + (row + 0.5) * gridCellSizeY;
                        
                        // 对于该 corner 的每个 channel，填充 MPV 值
                        for (size_t chIdx = 0; chIdx < cornerInfo.channels.size(); ++chIdx) {
                            int ch = cornerInfo.channels[chIdx];
                            Double_t mpv = mpvValues[padIdx][ch];
                            if (mpv > 0) {  // 只填充有效的 MPV 值
                                TH2D *th2d = cornerInfo.th2dChannels[chIdx];
                                th2d->Fill(padX, padY, mpv);
                            }
                        }
                    }
                }
            }
            
            mpvCornerCanvases.push_back(std::move(cornerInfo));
        }
        
        // 绘制 MPV TH2D 并标注数值
        for (auto &cornerInfo : mpvCornerCanvases) {
            for (size_t chIdx = 0; chIdx < cornerInfo.channels.size(); ++chIdx) {
                int ch = cornerInfo.channels[chIdx];
                cornerInfo.canvas->cd(chIdx + 1);
                gPad->SetTicks();
                TH2D *th2d = cornerInfo.th2dChannels[chIdx];
                th2d->Draw("COLZ");
                
                // 在 TH2D 上标注 MPV 值
                for (int binX = 1; binX <= cornerBlock; ++binX) {
                    for (int binY = 1; binY <= cornerBlock; ++binY) {
                        Double_t mpv = th2d->GetBinContent(binX, binY);
                        if (mpv > 0) {
                            // 获取 bin 的中心坐标（物理位置）
                            Double_t x = th2d->GetXaxis()->GetBinCenter(binX);
                            Double_t y = th2d->GetYaxis()->GetBinCenter(binY);
                            // 创建文本标签
                            TText *text = new TText(x, y, Form("%.0f", mpv));
                            text->SetTextAlign(22);  // 居中
                            text->SetTextSize(0.03);
                            text->SetTextColor(kBlack);
                            text->Draw();
                        }
                    }
                }
            }
        }
        
        TCanvas *canvas = new TCanvas("CellRealADCGrid", "CellRealADC grid", 2400, 2400);
        canvas->Divide(gridDivision, gridDivision);
        TDirectory *allPadDir = outputFile->mkdir("allPad");
        if (allPadDir) {
            allPadDir->cd();
        }
        struct CornerCanvasInfo {
            std::string name;
            TCanvas *canvas = nullptr;
            std::vector<int> padIndices;
        };
        std::vector<CornerCanvasInfo> cornerCanvases;
        std::map<int, std::pair<TCanvas*, int>> cornerPadSlots;
        std::map<int, TLegend*> cornerPadLegends;
        std::vector<TLegend*> cornerLegendsOwned;
        auto createCornerCanvas = [&](const std::string &name, int rowStart, int colStart) {
            CornerCanvasInfo info;
            info.name = name;
            info.canvas = new TCanvas(name.c_str(), name.c_str(), 1600, 1600);
            info.canvas->Divide(cornerBlock, cornerBlock);
            for (int dr = 0; dr < cornerBlock; ++dr) {
                for (int dc = 0; dc < cornerBlock; ++dc) {
                    int row = rowStart + dr;
                    int col = colStart + dc;
                    if (row < gridDivision && col < gridDivision) {
                        int padIdx = row * gridDivision + col;
                        info.padIndices.push_back(padIdx);
                        TLegend *cornerLeg = new TLegend(0.35, 0.45, 0.77, 0.75);
                        cornerLeg->SetBorderSize(0);
                        cornerLeg->SetFillStyle(0);
                        cornerLeg->SetTextSize(0.05);
                        cornerPadLegends[padIdx] = cornerLeg;
                        cornerLegendsOwned.push_back(cornerLeg);
                    }
                }
            }
            cornerCanvases.push_back(std::move(info));
        };
        createCornerCanvas("PadCorner_TopLeft", 0, 0);
        createCornerCanvas("PadCorner_TopRight", 0, gridDivision - cornerBlock);
        createCornerCanvas("PadCorner_BottomLeft", gridDivision - cornerBlock, 0);
        createCornerCanvas("PadCorner_BottomRight", gridDivision - cornerBlock,
                           gridDivision - cornerBlock);
        for (auto &corner : cornerCanvases) {
            for (size_t idx = 0; idx < corner.padIndices.size(); ++idx) {
                cornerPadSlots[corner.padIndices[idx]] =
                    {corner.canvas, static_cast<int>(idx)};
            }
        }
        std::vector<char> cornerPadDrawn(gridDivision * gridDivision, 0);
        for (int padIdx = 0; padIdx < gridDivision * gridDivision; ++padIdx) {
            canvas->cd(padIdx + 1);
            gPad->SetTicks();
            bool firstDrawn = false;
            bool firstDrawnPad = false;
            TCanvas *padCanvas = nullptr;
            TLegend *gridLegend = nullptr;
            TLegend *padLegend = nullptr;
            if (allPadDir) {
                std::string padCanvasName = "pad_" + std::to_string(padIdx);
                std::string padCanvasTitle = "Pad " + std::to_string(padIdx) + " overlays";
                padCanvas =
                    new TCanvas(padCanvasName.c_str(), padCanvasTitle.c_str(), 800, 800);
                padCanvas->cd();
                gPad->SetTicks();
                padLegend = new TLegend(0.58, 0.62, 0.9, 0.9);
                padLegend->SetBorderSize(0);
                padLegend->SetFillStyle(0);
                padLegend->SetTextSize(0.03);
            }
            for (int ch = 0; ch < 8; ++ch) {
                TH1D *hist = gridCellHistograms[padIdx][ch];
                if (!hist)
                    continue;
                hist->SetLineColor(colorPalette[ch % colorPalette.size()]);
                hist->SetLineWidth(2);
                hist->SetStats(0);
                hist->GetXaxis()->SetRangeUser(200., 3000.);
                hist->GetXaxis()->SetTitle("CellRealADC");
                hist->GetYaxis()->SetRangeUser(0., 50.);
                hist->GetYaxis()->SetTitle("Counts");
                hist->Draw(firstDrawn ? "HIST SAME" : "HIST");
                // 绘制拟合曲线
                TF1 *fitFunc = fitFunctions[padIdx][ch];
                if (fitFunc) {
                    fitFunc->Draw("SAME");
                }
                firstDrawn = true;
                if (!gridLegend) {
                    gridLegend = new TLegend(0.58, 0.62, 0.9, 0.9);
                    gridLegend->SetBorderSize(0);
                    gridLegend->SetFillStyle(0);
                    gridLegend->SetTextSize(0.03);
                }
                // 在 legend 中显示 MPV 值
                Double_t mpv = mpvValues[padIdx][ch];
                if (mpv > 0) {
                    gridLegend->AddEntry(hist, Form("ch%d MPV=%.0f", ch, mpv), "L");
                } else {
                    gridLegend->AddEntry(hist, Form("ch%d", ch), "L");
                }
                if (padCanvas) {
                    padCanvas->cd();
                    hist->Draw(firstDrawnPad ? "HIST SAME" : "HIST");
                    // 在 padCanvas 上也绘制拟合曲线
                    if (fitFunc) {
                        fitFunc->Draw("SAME");
                    }
                    if (padLegend) {
                        // 在 legend 中显示 MPV 值
                        Double_t mpv = mpvValues[padIdx][ch];
                        if (mpv > 0) {
                            padLegend->AddEntry(hist, Form("ch%d MPV=%.0f", ch, mpv), "L");
                        } else {
                            padLegend->AddEntry(hist, Form("ch%d", ch), "L");
                        }
                    }
                    firstDrawnPad = true;
                    canvas->cd(padIdx + 1);
                }
                auto cornerIt = cornerPadSlots.find(padIdx);
                if (cornerIt != cornerPadSlots.end()) {
                    TCanvas *cornerCanvas = cornerIt->second.first;
                    int slot = cornerIt->second.second;
                    if (cornerCanvas) {
                        cornerCanvas->cd(slot + 1);
                        gPad->SetTicks();
                        char &drawnFlag = cornerPadDrawn[padIdx];
                        hist->Draw(drawnFlag ? "HIST SAME" : "HIST");
                        // 在 cornerCanvas 上也绘制拟合曲线
                        if (fitFunc) {
                            fitFunc->Draw("SAME");
                        }
                        auto legIt = cornerPadLegends.find(padIdx);
                        if (legIt != cornerPadLegends.end() && legIt->second) {
                            // 在 legend 中显示 MPV 值
                            Double_t mpv = mpvValues[padIdx][ch];
                            if (mpv > 0) {
                                legIt->second->AddEntry(hist, Form("ch%d MPV=%.0f", ch, mpv), "L");
                            } else {
                                legIt->second->AddEntry(hist, Form("ch%d", ch), "L");
                            }
                        }
                        drawnFlag = 1;
                        canvas->cd(padIdx + 1);
                    }
                }
            }
            if (gridLegend && firstDrawn) {
                gridLegend->Draw();
            }
            if (padLegend && firstDrawnPad) {
                padCanvas->cd();
                padLegend->Draw();
                canvas->cd(padIdx + 1);
            }
            auto cornerLegendIt = cornerPadLegends.find(padIdx);
            auto cornerSlotIt = cornerPadSlots.find(padIdx);
            if (cornerLegendIt != cornerPadLegends.end() && cornerSlotIt != cornerPadSlots.end() &&
                cornerLegendIt->second && cornerPadDrawn[padIdx]) {
                TCanvas *cornerCanvas = cornerSlotIt->second.first;
                int slot = cornerSlotIt->second.second;
                cornerCanvas->cd(slot + 1);
                cornerLegendIt->second->Draw();
                canvas->cd(padIdx + 1);
            }
            if (padCanvas) {
                padCanvas->Write();
                delete padCanvas;
                if (allPadDir) {
                    allPadDir->cd();
                }
            }
        }
        canvas->SaveAs((outputDir + "/CellRealADC_grid.pdf").c_str());
        for (const auto &corner : cornerCanvases) {
            if (corner.canvas) {
                corner.canvas->SaveAs((outputDir + "/" + corner.name + ".pdf").c_str());
            }
        }
        // 保存 MPV canvas 为 PDF
        for (const auto &mpvCorner : mpvCornerCanvases) {
            if (mpvCorner.canvas) {
                mpvCorner.canvas->SaveAs((outputDir + "/" + mpvCorner.name + ".pdf").c_str());
            }
        }
        
        outputFile->cd();
        TDirectory *cornerDir = outputFile->mkdir("PadCornerCanvases");
        if (cornerDir) {
            cornerDir->cd();
            for (const auto &corner : cornerCanvases) {
                if (corner.canvas) {
                    corner.canvas->Write();
                }
            }
            outputFile->cd();
        }
        // 保存 MPV canvas 和 TH2D
        TDirectory *mpvDir = outputFile->mkdir("MPVCornerCanvases");
        if (mpvDir) {
            mpvDir->cd();
            for (const auto &mpvCorner : mpvCornerCanvases) {
                if (mpvCorner.canvas) {
                    mpvCorner.canvas->Write();
                }
                // 保存每个 channel 的 TH2D
                for (TH2D *th2d : mpvCorner.th2dChannels) {
                    if (th2d) {
                        th2d->Write();
                    }
                }
            }
            outputFile->cd();
        }
		TDirectory *histDirectory = outputFile->mkdir("CellRealADCHistograms");
        if (histDirectory) {
            histDirectory->cd();
			for (size_t padIdx = 0; padIdx < gridCellHistograms.size(); ++padIdx) {
				for (size_t ch = 0; ch < gridCellHistograms[padIdx].size(); ++ch) {
					TH1D *hist = gridCellHistograms[padIdx][ch];
					if (!hist)
						continue;
					hist->Write();
					// 保存拟合函数
					TF1 *fitFunc = fitFunctions[padIdx][ch];
					if (fitFunc) {
						fitFunc->Write();
					}
					delete hist;
				}
			}
        }
        if (xyDistribution) {
            outputFile->cd();
            xyDistribution->Write();
        }
        outputFile->cd();
        delete canvas;
        for (const auto &corner : cornerCanvases) {
            if (corner.canvas) {
                delete corner.canvas;
            }
        }
        for (TLegend *leg : cornerLegendsOwned) {
            delete leg;
        }
        // 清理 MPV canvas 和 TH2D
        for (auto &mpvCorner : mpvCornerCanvases) {
            if (mpvCorner.canvas) {
                delete mpvCorner.canvas;
            }
            for (TH2D *th2d : mpvCorner.th2dChannels) {
                if (th2d) {
                    delete th2d;
                }
            }
        }
        delete xyDistribution;
        outputFile->Close();
        delete outputFile;
        // 清理拟合函数
        for (auto &padVec : fitFunctions) {
            for (TF1 *fitFunc : padVec) {
                if (fitFunc) {
                    delete fitFunc;
                }
            }
        }
        std::cout << "匹配数据已保存至 " << outputFileName << "，共写入 " << savedEntries << " 条记录。" << std::endl;
    } else {
        std::cout << "未找到需要保存的数据，未创建输出文件。" << std::endl;
    }
    
}
