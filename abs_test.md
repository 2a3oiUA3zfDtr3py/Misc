#### 絶対値を求める方法の実行時間を調べた
環境  
Linux localhost 4.10.6-gentoo #1 SMP　x86_64 Intel(R) Core(TM)2 Duo CPU P8600 @ 2.40GHz GenuineIntel GNU/Linux  
gcc version 5.4.0  
ソース  
<https://github.com/2a3oiUA3zfDtr3py/Misc/blob/master/abs_test.c>  
コンパイルオプション  
-O3 -D type?  

三回実行した時間の中央値
type1 user    0m10.910s  
type2 user    0m10.940s  
type3 user    0m4.560s  
type4 user    0m11.000s  
type5 user    0m10.910s  
type6 user    0m4.570s  
type7 user    0m10.990s  
type8 user    0m4.620s  

andで符号ビットを0で上書きするだけのtype1が一番早そうなのに、  
if文の条件式の中身がビットシフトか0との比較の時にしか早くならない。  
type3とtype6はマシだが普通に書いたtype8と差が殆ど無い。  
&gt;&gt;は早いが&amp;は遅いようだ。  
