# 2511_mm_csi / 简介

## micromegas_ana Quick Guide / 快速指南
```
cd micromegas_ana
./run.sh
```
- Inputs / 输入：`mmdata/`, `decodedata/` → 合并成 `result/adas_track_data.root`、`result/result_decode.root`。
- Outputs / 输出：全部写入 `result/`（含 `common_trigger_data.root`, `CellRealADC_grid.pdf`, `PadCorner_*.pdf`, per-pad canvases）。
- Needs ROOT 6.x (`root`, `hadd`). 运行前请确认已安装 ROOT 6.x。

## vdecode Quick Guide / 快速指南
```
cd vdecode
mkdir -p build 
cd build && cmake ..
make -j12 && cd ..
source run.sh        # or run decoder manually / 也可自行执行解码程序
```
- Sources & configs: `src/`, `include/`; raw data in `data/`。结果写入 `vdecode/result/` 并为 `micromegas_ana` 所用。

## Notes / 说明
- Keep large intermediates inside each `result/` directory; scripts assume relative paths. 大文件统一放在各自 `result/` 目录。
- Update this README when workflows/datasets change. 如流程或数据集有更新，请同步修改 README。

