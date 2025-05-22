# petagenome pipeline

## ダウンロード方法

```
$ git clone https://github.com/jinn0027/petagenome_pipeline.git
```

## ホスト環境設定

以下のスクリプトはalmalinux9にて実行確認。ubuntu等では修正が必要。

```
$ cd petagenome_pipeline/etc
$ source ./host_setup.sh
```

ここでは以下を行っている。

- apptainerインストール
- nextflowインストール
- sifファイルからsandboxを経由せずに直接実行するためにsquashfuseをインストール
- apptainerのfakerootでdnfを行うためにはselinuxをpermissiveに設定

## 外部リソースダウンロード

以下のスクリプトはalmalinux9にて実行確認。ubuntu等では修正が必要。
外部のコードやデータ等をダウンロードする。

```
$ cd petagenome_pipeline/etc
$ ./download.sh
```

ダウンロード先は petagenome_pipeline/external となる。

## モジュール

モジュール名 foo の場合、モジュールは petagenome_pipeline/modules/foo 以下にある。
モジュール一覧は[こちら](doc/modules.md)。

最初に各モジュールが利用する共有のsifファイルを作成する必要がある。

```
$ cd petagenome_pipeline/modules/common
$ make
```

次に以下の手順でsandboxとコンテナをfoo.sbx、foo.sif という名前で作成する。

```
$ cd petagenome_pipeline/modules/foo
$ make
```

モジュール毎の各ディレクトリにはtestサブディレクトリがあり、ここで簡易テストを行える。

```
$ cd petagenome_pipeline/modules/foo/test
$ ./t.sh
```



