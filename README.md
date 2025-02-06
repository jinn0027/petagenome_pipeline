# petagenome pipeline

## ダウンロード方法

`
$ git clone https://github.com/jinn0027/petagenome_pipeline.git
`

## ホスト環境設定

almalinux9にて実行。

```$ cd petagenome_pipeline/etc
$ ./host_setup.sh
```

ここでは以下を行っている。

- apptainerインストール
- nextflowインストール
- sifファイルからsandboxを経由せずに直接実行するためにsquashfuseをインストール
- apptainerのfakerootでdnfを行うためにはselinuxをpermissiveに設定

## 外部リソースダウンロード

almalinux9にて実行。

```$ cd petagenome_pipeline/etc
$ ./download.sh
```

## モジュールテスト




