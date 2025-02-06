# petagenome pipeline

## ダウンロード方法

  $ git clone https://github.com/jinn0027/petagenome_pipeline.git

## ホスト環境設定

  almalinux9にて実行。

  $ cd petagenome_pipeline/etc
  $ ./host_setup.sh

  ここでは以下を行っている。

  ### nextflowとapptainerの設定を行っている。
  ### またapptainerのfakerootでdnfを行うためにはSELinuxをpermissiveにしている。
  ### sifファイルから実行するためにsquashfuseをインストール

## 外部リソースダウンロード

  almalinux9にて実行。

  $ cd petagenome_pipeline/etc
  $ ./download.sh

## モジュールテスト




