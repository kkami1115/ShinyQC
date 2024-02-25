# ShinyQC
指定したディレクトリ内のfastqファイルを自動でペアにして、Rfastpをpaired-endで実行するGUIです。

# 使い方
`fastq_dir` ... fastqが存在するディレクトリ。

`result_dir` ... 結果(html,json)と処理済fastqを出力するディレクトリ。

`file_pattern` ... Read1とRead2のパターンを正規表現で。例: `"_R[12]"`

全ファイルの処理が終了するとOverviewタブとVisualize per sampleタブの内容が表示できるようになります。


# ToDo
- `furrr`の非同期処理で立てているプロセスのコンソール出力をShinyにログとして出力する。

