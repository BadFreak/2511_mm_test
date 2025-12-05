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
#include "TLatex.h"

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
        std::vector<Double_t> outCellRealADC;
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
                    // 计算角度修正因子：sqrt(1 + x_slope² + y_slope²)
                    Double_t lengthFactor = TMath::Sqrt(1.0 + outXSlope * outXSlope + outYSlope * outYSlope);
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
                        // 计算原始 realADC
                        Double_t realADC = static_cast<Double_t>(outCellADC[idx] - outCellPLAT[idx]);
                        // 应用角度修正：除以路径长度因子
                        realADC /= lengthFactor;
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
        
        // 创建8个TH2D，每个channel对应其corner区域的MPV分布
        std::vector<TH2D*> mpvChannelTH2D(8, nullptr);
        
        // 定义每个channel的corner信息
        struct ChannelCornerInfo {
            int rowStart;
            int colStart;
            Double_t xMin, xMax;
            Double_t yMin, yMax;
            std::string cornerName;
        };
        
        std::vector<ChannelCornerInfo> channelInfo(8);
        // TopRight (col=0-4, row=0-4): channels 0, 1
        channelInfo[0] = {0, 0, gridXMin, gridXMin + cornerBlock * gridCellSizeX, 
                         gridYMin, gridYMin + cornerBlock * gridCellSizeY, "TopRight"};
        channelInfo[1] = channelInfo[0];
        
        // TopLeft (col=5-9, row=0-4): channels 2, 3
        channelInfo[2] = {0, cornerBlock, gridXMin + cornerBlock * gridCellSizeX, gridXMax,
                         gridYMin, gridYMin + cornerBlock * gridCellSizeY, "TopLeft"};
        channelInfo[3] = channelInfo[2];
        
        // BottomRight (col=0-4, row=5-9): channels 4, 5
        channelInfo[4] = {cornerBlock, 0, gridXMin, gridXMin + cornerBlock * gridCellSizeX,
                         gridYMin + cornerBlock * gridCellSizeY, gridYMax, "BottomRight"};
        channelInfo[5] = channelInfo[4];
        
        // BottomLeft (col=5-9, row=5-9): channels 6, 7
        channelInfo[6] = {cornerBlock, cornerBlock, gridXMin + cornerBlock * gridCellSizeX, gridXMax,
                         gridYMin + cornerBlock * gridCellSizeY, gridYMax, "BottomLeft"};
        channelInfo[7] = channelInfo[6];
        
        // 创建并填充TH2D
        for (int ch = 0; ch < 8; ++ch) {
            std::string histName = "MPV_Channel" + std::to_string(ch);
            std::string histTitle = "Channel " + std::to_string(ch) + " (" + channelInfo[ch].cornerName + 
                                   ") MPV;X position (mm);Y position (mm);MPV";
            
            // 创建5×5的bin，对应corner的25个pad
            mpvChannelTH2D[ch] = new TH2D(histName.c_str(), histTitle.c_str(),
                                          cornerBlock, channelInfo[ch].xMin, channelInfo[ch].xMax,
                                          cornerBlock, channelInfo[ch].yMin, channelInfo[ch].yMax);
            mpvChannelTH2D[ch]->SetStats(0);
            
            // 遍历该corner的25个pad，填充MPV值
            int rowStart = channelInfo[ch].rowStart;
            int colStart = channelInfo[ch].colStart;
            
            for (int dr = 0; dr < cornerBlock; ++dr) {
                for (int dc = 0; dc < cornerBlock; ++dc) {
                    int row = rowStart + dr;
                    int col = colStart + dc;
                    
                    // 计算padIdx
                    int colForCanvas = (gridDivision - 1) - col;
                    int padIdx = row * gridDivision + colForCanvas;
                    
                    // 获取MPV值
                    Double_t mpv = mpvValues[padIdx][ch];
                    
                    if (mpv > 0) {  // 只填充有效的MPV值
                        // 计算该pad的中心坐标
                        Double_t padX = gridXMin + (col + 0.5) * gridCellSizeX;
                        Double_t padY = gridYMin + (row + 0.5) * gridCellSizeY;
                        
                        // 填充到TH2D
                        mpvChannelTH2D[ch]->Fill(padX, padY, mpv);
                    }
                }
            }
        }
        
        // 创建8页PDF，每页显示一个channel的MPV分布
        TCanvas *mpvChannelCanvas = new TCanvas("MPVChannelCanvas", "MPV by Channel", 1200, 1000);
        std::string pdfFileName = outputDir + "/MPV_by_Channel.pdf";
        
        for (int ch = 0; ch < 8; ++ch) {
            mpvChannelCanvas->cd();
            mpvChannelCanvas->Clear();
            gPad->SetTicks();
            gPad->SetRightMargin(0.18);  // 增大右边距，为color bar预留空间
            gPad->SetLeftMargin(0.12);
            gPad->SetTopMargin(0.12);    // 增大上边距，避免title和color bar重叠
            gPad->SetBottomMargin(0.12);
            
            if (mpvChannelTH2D[ch]) {
                // 设置bin内容显示格式为整数（无小数）
                mpvChannelTH2D[ch]->SetMarkerSize(1.8);  // 增大TEXT字体大小
                gStyle->SetPaintTextFormat("1.0f");     // 设置显示格式为整数
                
                mpvChannelTH2D[ch]->Draw("COLZ TEXT");  // TEXT选项显示每个bin的数值
                
                // 调整Z轴title位置，避免与color bar的label重叠
                mpvChannelTH2D[ch]->GetZaxis()->SetTitleOffset(1.8);  // 增大offset
                mpvChannelTH2D[ch]->GetZaxis()->SetTitleSize(0.04);    // 稍微减小title字体
            }
            
            // 保存PDF
            if (ch == 0) {
                mpvChannelCanvas->Print((pdfFileName + "(").c_str());
            } else if (ch == 7) {
                mpvChannelCanvas->Print((pdfFileName + ")").c_str());
            } else {
                mpvChannelCanvas->Print(pdfFileName.c_str());
            }
        }
        
        std::cout << "\nMPV by Channel PDF 已保存至: " << pdfFileName << std::endl;
        
        // 创建一个8 pad的canvas，按照corner位置排列
        TCanvas *mpvAllChannelsCanvas = new TCanvas("MPVAllChannels", "MPV All Channels", 2400, 1200);
        mpvAllChannelsCanvas->Divide(4, 2, 0.01, 0.01);  // 4列2行，小间隙
        
        // 定义8个pad的绘图顺序和对应的channel
        // 布局：
        // 1(ch2)  2(ch3)  |  3(ch0)  4(ch1)     <- TopLeft  | TopRight
        // 5(ch6)  6(ch7)  |  7(ch4)  8(ch5)     <- BottomLeft | BottomRight
        std::vector<int> channelOrder = {2, 3, 0, 1, 6, 7, 4, 5};
        
        for (int i = 0; i < 8; ++i) {
            int ch = channelOrder[i];
            mpvAllChannelsCanvas->cd(i + 1);
            gPad->SetTicks();
            gPad->SetRightMargin(0.15);
            gPad->SetLeftMargin(0.12);
            gPad->SetTopMargin(0.12);
            gPad->SetBottomMargin(0.12);
            
            if (mpvChannelTH2D[ch]) {
                mpvChannelTH2D[ch]->SetMarkerSize(1.8);
                gStyle->SetPaintTextFormat("1.0f");
                mpvChannelTH2D[ch]->Draw("COLZ TEXT");
                mpvChannelTH2D[ch]->GetZaxis()->SetTitleOffset(1.8);  // 增大offset
                mpvChannelTH2D[ch]->GetZaxis()->SetTitleSize(0.04);    // 稍微减小title字体
            }
        }
        
        // 保存为单页PDF
        std::string allChannelsPdfFileName = outputDir + "/MPV_All_Channels.pdf";
        mpvAllChannelsCanvas->SaveAs(allChannelsPdfFileName.c_str());
        std::cout << "MPV All Channels PDF 已保存至: " << allChannelsPdfFileName << std::endl;
        
        // 创建8个TH1D，收集每个channel的MPV分布
        std::vector<TH1D*> mpvChannelTH1D(8, nullptr);
        
        for (int ch = 0; ch < 8; ++ch) {
            std::string histName = "MPV_Hist_Channel" + std::to_string(ch);
            std::string histTitle = "Channel " + std::to_string(ch) + " (" + channelInfo[ch].cornerName + 
                                   ") MPV Distribution;MPV;Entries";
            
            // 创建TH1D，范围根据MPV的典型值设置
            mpvChannelTH1D[ch] = new TH1D(histName.c_str(), histTitle.c_str(), 50, 700., 1900.);
            mpvChannelTH1D[ch]->SetStats(0);  // 不显示stat box
            mpvChannelTH1D[ch]->SetLineColor(kBlack);  // 黑线
            mpvChannelTH1D[ch]->SetLineWidth(2);
            mpvChannelTH1D[ch]->SetFillColor(kBlue);   // 蓝色填充
            mpvChannelTH1D[ch]->SetFillStyle(3004);    // 花纹样式
            
            // 从TH2D中提取所有有效的MPV值填充到TH1D
            if (mpvChannelTH2D[ch]) {
                for (int binX = 1; binX <= cornerBlock; ++binX) {
                    for (int binY = 1; binY <= cornerBlock; ++binY) {
                        Double_t mpv = mpvChannelTH2D[ch]->GetBinContent(binX, binY);
                        if (mpv > 0) {  // 只填充有效的MPV值
                            mpvChannelTH1D[ch]->Fill(mpv);
                        }
                    }
                }
            }
        }
        
        // 创建8 pad canvas显示MPV分布histogram
        TCanvas *mpvHistCanvas = new TCanvas("MPVHistCanvas", "MPV Distribution by Channel", 2400, 1200);
        mpvHistCanvas->Divide(4, 2, 0.01, 0.01);  // 4列2行，和TH2D的布局一致
        
        for (int i = 0; i < 8; ++i) {
            int ch = channelOrder[i];
            mpvHistCanvas->cd(i + 1);
            gPad->SetTicks();
            gPad->SetLeftMargin(0.12);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.1);
            gPad->SetBottomMargin(0.12);
            
            if (mpvChannelTH1D[ch]) {
                mpvChannelTH1D[ch]->Draw("HIST");
                
                // 计算统计量
                if (mpvChannelTH1D[ch]->GetEntries() > 0) {
                    Double_t mean = mpvChannelTH1D[ch]->GetMean();
                    Double_t rms = mpvChannelTH1D[ch]->GetRMS();
                    Double_t rmsOverMean = (mean > 0) ? rms / mean : 0.0;
                    
                    // 使用TLatex显示统计信息
                    TLatex *latex = new TLatex();
                    latex->SetNDC(true);
                    latex->SetTextFont(62);         // 粗体字体，更明显
                    latex->SetTextSize(0.06);       // 增大字体
                    latex->SetTextColor(kRed + 1);  // 红色，更醒目
                    
                    // 在图上显示MEAN, RMS, RMS/MEAN（向左下移动）
                    Double_t textX = 0.45;   // 从0.55向左移到0.45
                    Double_t textY1 = 0.75;  // 从0.85向下移到0.75
                    Double_t textY2 = 0.65;  // 从0.75向下移到0.65
                    Double_t textY3 = 0.55;  // 从0.65向下移到0.55
                    
                    latex->DrawLatex(textX, textY1, Form("Mean = %.1f", mean));
                    latex->DrawLatex(textX, textY2, Form("RMS = %.1f", rms));
                    latex->DrawLatex(textX, textY3, Form("RMS/Mean = %.4f", rmsOverMean));
                }
            }
        }
        
        // 保存为单页PDF
        std::string mpvHistPdfFileName = outputDir + "/MPV_Distribution_All_Channels.pdf";
        mpvHistCanvas->SaveAs(mpvHistPdfFileName.c_str());
        std::cout << "MPV Distribution All Channels PDF 已保存至: " << mpvHistPdfFileName << std::endl;
        
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
            // 收集该corner的所有padIdx，然后按照padIdx顺序排列
            // 这样PadCorner的排列和CellRealADC_grid完全一致
            std::vector<int> cornerPadIndices;
            for (int dr = 0; dr < cornerBlock; ++dr) {
                for (int dc = 0; dc < cornerBlock; ++dc) {
                    int row = rowStart + dr;
                    int col = colStart + dc;
                    if (row < gridDivision && col < gridDivision) {
                        int colForCanvas = (gridDivision - 1) - col;
                        int padIdx = row * gridDivision + colForCanvas;
                        cornerPadIndices.push_back(padIdx);
                    }
                }
            }
            // 按padIdx从小到大排序，保持和CellRealADC_grid相同的顺序
            std::sort(cornerPadIndices.begin(), cornerPadIndices.end());
            
            for (int padIdx : cornerPadIndices) {
                info.padIndices.push_back(padIdx);
                TLegend *cornerLeg = new TLegend(0.35, 0.45, 0.77, 0.75);
                cornerLeg->SetBorderSize(0);
                cornerLeg->SetFillStyle(0);
                cornerLeg->SetTextSize(0.05);
                cornerPadLegends[padIdx] = cornerLeg;
                cornerLegendsOwned.push_back(cornerLeg);
            }
            cornerCanvases.push_back(std::move(info));
        };
        // 注意：原点在右上角，X轴正向向左，Y轴正向向下
        // col=0-4 对应 X=160-310 (右侧), col=5-9 对应 X=310-460 (左侧)
        // row=0-4 对应 Y=140.5-290.5 (上侧), row=5-9 对应 Y=290.5-440.5 (下侧)
        createCornerCanvas("PadCorner_TopRight", 0, 0);                              // 右上: col=0-4, row=0-4
        createCornerCanvas("PadCorner_TopLeft", 0, gridDivision - cornerBlock);      // 左上: col=5-9, row=0-4
        createCornerCanvas("PadCorner_BottomRight", gridDivision - cornerBlock, 0);  // 右下: col=0-4, row=5-9
        createCornerCanvas("PadCorner_BottomLeft", gridDivision - cornerBlock,       // 左下: col=5-9, row=5-9
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
                hist->GetYaxis()->SetRangeUser(0., 60.);
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
        // 保存8个channel的MPV TH2D和TH1D
        TDirectory *mpvChannelDir = outputFile->mkdir("MPVbyChannel");
        if (mpvChannelDir) {
            mpvChannelDir->cd();
            for (int ch = 0; ch < 8; ++ch) {
                if (mpvChannelTH2D[ch]) {
                    mpvChannelTH2D[ch]->Write();
                }
                if (mpvChannelTH1D[ch]) {
                    mpvChannelTH1D[ch]->Write();
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
        // 清理MPV channel canvas和TH2D/TH1D
        if (mpvChannelCanvas) {
            delete mpvChannelCanvas;
        }
        if (mpvAllChannelsCanvas) {
            delete mpvAllChannelsCanvas;
        }
        if (mpvHistCanvas) {
            delete mpvHistCanvas;
        }
        for (int ch = 0; ch < 8; ++ch) {
            if (mpvChannelTH2D[ch]) {
                delete mpvChannelTH2D[ch];
            }
            if (mpvChannelTH1D[ch]) {
                delete mpvChannelTH1D[ch];
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
