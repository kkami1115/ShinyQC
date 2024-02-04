# ShinyQC
指定したディレクトリ内のfastqをRfastpを使用してpairedで実行するGUIです。

# 使い方
`fastq_dir` ... fastqが存在するディレクトリ。

`result_dir` ... 結果(html,json)と処理済fastqを出力するディレクトリ。

`file_pattern` ... Read1とRead2のパターンを正規表現で。例: `"_R[12]"`

全ファイルの処理が終了するとOverviewタブとVisualize per sampleタブの内容が表示できるようになります。


# ToDo
- 実行〜処理終了の通知
- `furrr`の非同期処理で立てているプロセスのコンソール出力をShinyにログとして出力する。
- `base`/`ggplot2`の図を`plotly`化

