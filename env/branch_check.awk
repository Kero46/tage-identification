# branch_check.awk — 逆アセンブルから「データ依存の条件分岐」を数える
#
# 目的:
#   空 asm やコンパイラ属性が効いていることを前提にせず、実際に生成された
#   機械語で不変条件を検査する。
#     §1.2 条件分岐が消えていない（if-conversion されていない）
#     §1.3 論理的な 1 サイトが 1 つの PC に対応する（ブロック複製が無い）
#
# 判定の考え方:
#   パターン配列の 1 バイト読み出しを 1 サイトの起点とみなし、その値を
#   使う条件分岐が存在することをレジスタの一致で追跡する。
#   レジスタを突き合わせるので、ループ制御の分岐は数えない。ループは
#   ループカウンタのレジスタを使うため、パターン値のレジスタと一致しない。
#   （シンボル名だけ直して条件分岐の有無を grep すると、ループ制御の
#   b.eq / jne を拾って誤って通る。実際にその危険があった。）
#   追跡の途中に条件選択命令（cmov / csel / csinc 等）が現れたら、その
#   サイトは if-conversion されている。
#
# 入力: objdump -d の出力全体
# 変数: SYM=検査するシンボル名  ARCH=x86|arm64  EXPECT=期待サイト数
# 出力: 1 行
#   found=<0|1> srcs=<読み出し数> sites=<データ依存分岐> sel=<条件選択>
#   cond=<条件分岐総数> ng=<理由|->
#
# 検査不能（シンボルが無い・命令を 1 つも解釈できない）は成功にしない。
# 呼び出し側が found=0 や srcs=0 を失敗として扱う。

# ---------------------------------------------------------------- 字句
function ishex(s,   i, c) {
    if (s == "") return 0
    for (i = 1; i <= length(s); i++) {
        c = substr(s, i, 1)
        if (index("0123456789abcdefABCDEF", c) == 0) return 0
    }
    return 1
}
# 命令バイト列の桁: GNU objdump は 2 桁区切り、llvm-objdump の arm64 は 8 桁。
# 語長 3 の "add" 等と衝突しないよう長さを 2 と 8 に限る。
function isbytecol(s) { return ishex(s) && (length(s) == 2 || length(s) == 8) }

function trim(s) { gsub(/^[ \t]+|[ \t]+$/, "", s); return s }

# ---------------------------------------------------------------- レジスタ正規化
# x86: %rax/%eax/%ax/%al/%ah を同一視する（バイト別名で値を使うため必須）
function canon_x86(r) {
    if (r ~ /^r[0-9]+[dwb]?$/) { sub(/[dwb]$/, "", r); return r }
    if (r ~ /^(rax|eax|ax|al|ah)$/) return "a"
    if (r ~ /^(rbx|ebx|bx|bl|bh)$/) return "b"
    if (r ~ /^(rcx|ecx|cx|cl|ch)$/) return "c"
    if (r ~ /^(rdx|edx|dx|dl|dh)$/) return "d"
    if (r ~ /^(rsi|esi|si|sil)$/)   return "si"
    if (r ~ /^(rdi|edi|di|dil)$/)   return "di"
    if (r ~ /^(rbp|ebp|bp|bpl)$/)   return "bp"
    if (r ~ /^(rsp|esp|sp|spl)$/)   return "sp"
    return ""
}
# arm64: w12 と x12 は同一レジスタ
function canon_arm(r) {
    if (r ~ /^[wx][0-9]+$/) return substr(r, 2)
    if (r ~ /^[wx]zr$/)     return "zr"
    if (r == "sp")          return "sp"
    return ""
}
function canon(r) { return (ARCH == "x86") ? canon_x86(r) : canon_arm(r) }

# 命令が言及する全レジスタを "|a|d|" 形式で返す
function regids(opstr,   arr, n, i, out, r, t) {
    out = "|"
    n = split(opstr, arr, /[^A-Za-z0-9%]+/)
    for (i = 1; i <= n; i++) {
        t = arr[i]
        if (ARCH == "x86") {
            if (t !~ /^%/) continue
            r = canon(substr(t, 2))
        } else {
            r = canon(t)
        }
        if (r != "") out = out r "|"
    }
    return out
}
function mentions(k, reg) { return index(regids(ops[k]), "|" reg "|") > 0 }

# 第 1 オペランド（arm64 の宛先）
function first_op(k,   arr, n) {
    n = split(ops[k], arr, ",")
    return canon(trim(arr[1]))
}
# 最終オペランド（x86 AT&T の宛先）
function last_op(k,   arr, n, t) {
    n = split(ops[k], arr, ",")
    t = trim(arr[n])
    sub(/^%/, "", t)
    return canon(t)
}

# ---------------------------------------------------------------- 命令分類
function is_select(k) {
    if (ARCH == "x86") return m[k] ~ /^cmov/
    # csel/csinc/csinv/csneg と別名（cset/cinc/…）。ccmp/ccmn は条件を
    # 1 本の分岐に畳み込むため、これが出たらサイトが融合している。
    return m[k] ~ /^(csel|csinc|csinv|csneg|cset|csetm|cinc|cinv|cneg|ccmp|ccmn)$/
}
function is_condbr(k) {
    if (ARCH == "x86") return m[k] ~ /^j(e|ne|z|nz|s|ns|a|ae|b|be|g|ge|l|le|p|np|o|no|c|nc)$/
    return m[k] ~ /^b\./ || m[k] ~ /^(cbz|cbnz|tbz|tbnz)$/
}
# フラグを立てる比較（第 1 オペランドが被検査レジスタか見る）
function is_cmp(k) {
    if (ARCH == "x86") return m[k] ~ /^(cmp|test|cmpb|testb|cmpl|testl|cmpq|testq|cmpw|testw)$/
    return m[k] ~ /^(cmp|cmn|tst|subs|adds|ands|bics|negs)$/
}
# 第 1 オペランドを書き換えるか（arm64）
function writes_first(k) {
    return m[k] !~ /^b\./ &&
           m[k] !~ /^(str|strb|strh|stur|sturb|sturh|stp|stlr|stlrb|cmp|cmn|tst|ccmp|ccmn|ret|nop|b|bl|br|blr|cbz|cbnz|tbz|tbnz|dmb|dsb|isb|brk|svc)$/
}
# パターン配列の 1 バイト読み出しか
function is_byte_src(k) {
    if (ARCH == "x86") {
        # メモリ参照を伴うバイト読み出し、またはメモリを直接見るバイト比較
        if (ops[k] !~ /\(/) return 0
        return m[k] ~ /^(movzbl|movzbq|movzbw|movsbl|movsbq|movsbw|movb)$/ ||
               m[k] ~ /^(cmpb|testb)$/
    }
    return m[k] ~ /^(ldrb|ldurb|ldrsb|ldursb)$/
}
# 比較そのものが読み出しを兼ねる形（cmpb $0x0,(%rax)）か
function src_sets_flags(k) { return ARCH == "x86" && m[k] ~ /^(cmpb|testb)$/ }

# ---------------------------------------------------------------- 本体の切り出し
/^[0-9a-fA-F]+[ \t]+<.*>:[ \t]*$/ {
    name = $0
    sub(/^[0-9a-fA-F]+[ \t]+</, "", name)
    sub(/>:[ \t]*$/, "", name)
    sub(/^_/, "", name)          # Mach-O のシンボル接頭辞
    inbody = (name == SYM)
    if (inbody) found = 1
    next
}
inbody && /^[ \t]*$/ { inbody = 0; next }
inbody {
    line = $0
    sub(/[;#].*$/, "", line)                             # 注釈を落とす
    sub(/^[ \t]*[0-9a-fA-F]+:[ \t]*/, "", line)          # 番地を落とす
    nf = split(line, f, /[ \t]+/)
    mn = ""; rest = ""
    for (i = 1; i <= nf; i++) {
        if (f[i] == "") continue
        if (mn == "") { if (isbytecol(f[i])) continue; mn = f[i]; continue }
        rest = (rest == "") ? f[i] : rest " " f[i]
    }
    if (mn == "") next
    n++
    m[n] = mn
    ops[n] = rest
    next
}

# ---------------------------------------------------------------- 判定
END {
    if (!found) { print "found=0 srcs=0 sites=0 sel=0 cond=0 ng=シンボルが見つからない"; exit }
    if (n == 0) { print "found=1 srcs=0 sites=0 sel=0 cond=0 ng=命令を解釈できない"; exit }

    WINDOW = 16          # 1 サイトは読み出しから分岐まで数命令。余裕を見て 16
    srcs = 0; sites = 0; sel = 0; cond = 0; broken = 0
    for (k = 1; k <= n; k++) {
        if (is_select(k)) sel++
        if (is_condbr(k)) cond++
    }
    for (k = 1; k <= n; k++) {
        if (!is_byte_src(k)) continue
        srcs++
        reg = (ARCH == "x86") ? last_op(k) : first_op(k)
        flagged = src_sets_flags(k)          # cmpb $0,(mem) は自分でフラグを立てる
        hit = 0
        for (j = k + 1; j <= n && j <= k + WINDOW; j++) {
            if (is_select(j)) { broken++; break }         # if-conversion された
            if (ARCH != "x86" && m[j] ~ /^(cbz|cbnz|tbz|tbnz)$/) {
                if (first_op(j) == reg) { hit = 1; break }
                continue
            }
            if (is_cmp(j) && (flagged == 0)) {
                if (ARCH == "x86") { if (reg != "" && mentions(j, reg)) flagged = 1 }
                else               { if (first_op(j) == reg)           flagged = 1 }
                continue
            }
            if (is_condbr(j)) {
                if (flagged) { hit = 1; break }
                continue                                   # 他の条件（ループ等）
            }
            # 値が使われる前に上書きされたら追跡打ち切り
            if (ARCH == "x86") { if (reg != "" && last_op(j) == reg && m[j] ~ /^(mov|lea|add|sub|xor|and|or)/) break }
            else               { if (writes_first(j) && first_op(j) == reg) break }
        }
        if (hit) sites++
    }

    ng = "-"
    if (sel > 0)            ng = sprintf("条件選択命令が %d 個（if-conversion）", sel)
    else if (broken > 0)    ng = sprintf("%d サイトが条件選択に潰れている", broken)
    else if (srcs == 0)     ng = "パターン読み出しを検出できない"
    else if (srcs != EXPECT) ng = sprintf("読み出し %d 個（期待 %d）", srcs, EXPECT)
    else if (sites != EXPECT) ng = sprintf("データ依存分岐 %d 個（期待 %d）", sites, EXPECT)
    printf "found=1 srcs=%d sites=%d sel=%d cond=%d ng=%s\n", srcs, sites, sel, cond, ng
}
