# ShinyQC
指定したディレクトリ内のfastqファイルを自動でペアにして、Rfastpをpaired-endで実行するGUIです。

# 使い方
`fastq_dir` ... fastqが存在するディレクトリ。

`result_dir` ... 結果(html,json)と処理済fastqを出力するディレクトリ。

`file_pattern` ... Read1とRead2のパターンを正規表現で。例: `"_R[12]"`

全ファイルの処理が終了するとOverviewタブとVisualize per sampleタブの内容が表示できるようになります。


# ToDo
- per sampleの可視化を実際のhtmlと合わせる（今は微妙に仕様が違う）
- duplicate rateをOverviewで確認できるようにする
        - `future`のplanをユーザー側が変更できるようにするとか、`Rfastp`に渡す引数をどうするのとか、設定タブを整備する