# Awk script for OpenFOAM log file extraction
BEGIN {
    Iteration=0
    resetCounters()
}

# Reset counters used for variable postfix
function resetCounters() {
    clockTimeCnt=0
    executionTimeCnt=0
    ILambda_0_0Cnt=0
    ILambda_0_0FinalResCnt=0
    ILambda_0_0ItersCnt=0
    ILambda_1_0Cnt=0
    ILambda_10_0Cnt=0
    ILambda_10_0FinalResCnt=0
    ILambda_10_0ItersCnt=0
    ILambda_1_0FinalResCnt=0
    ILambda_1_0ItersCnt=0
    ILambda_11_0Cnt=0
    ILambda_11_0FinalResCnt=0
    ILambda_11_0ItersCnt=0
    ILambda_12_0Cnt=0
    ILambda_12_0FinalResCnt=0
    ILambda_12_0ItersCnt=0
    ILambda_13_0Cnt=0
    ILambda_13_0FinalResCnt=0
    ILambda_13_0ItersCnt=0
    ILambda_14_0Cnt=0
    ILambda_14_0FinalResCnt=0
    ILambda_14_0ItersCnt=0
    ILambda_15_0Cnt=0
    ILambda_15_0FinalResCnt=0
    ILambda_15_0ItersCnt=0
    ILambda_16_0Cnt=0
    ILambda_16_0FinalResCnt=0
    ILambda_16_0ItersCnt=0
    ILambda_17_0Cnt=0
    ILambda_17_0FinalResCnt=0
    ILambda_17_0ItersCnt=0
    ILambda_18_0Cnt=0
    ILambda_18_0FinalResCnt=0
    ILambda_18_0ItersCnt=0
    ILambda_19_0Cnt=0
    ILambda_19_0FinalResCnt=0
    ILambda_19_0ItersCnt=0
    ILambda_2_0Cnt=0
    ILambda_20_0Cnt=0
    ILambda_20_0FinalResCnt=0
    ILambda_20_0ItersCnt=0
    ILambda_2_0FinalResCnt=0
    ILambda_2_0ItersCnt=0
    ILambda_21_0Cnt=0
    ILambda_21_0FinalResCnt=0
    ILambda_21_0ItersCnt=0
    ILambda_22_0Cnt=0
    ILambda_22_0FinalResCnt=0
    ILambda_22_0ItersCnt=0
    ILambda_23_0Cnt=0
    ILambda_23_0FinalResCnt=0
    ILambda_23_0ItersCnt=0
    ILambda_24_0Cnt=0
    ILambda_24_0FinalResCnt=0
    ILambda_24_0ItersCnt=0
    ILambda_25_0Cnt=0
    ILambda_25_0FinalResCnt=0
    ILambda_25_0ItersCnt=0
    ILambda_26_0Cnt=0
    ILambda_26_0FinalResCnt=0
    ILambda_26_0ItersCnt=0
    ILambda_27_0Cnt=0
    ILambda_27_0FinalResCnt=0
    ILambda_27_0ItersCnt=0
    ILambda_28_0Cnt=0
    ILambda_28_0FinalResCnt=0
    ILambda_28_0ItersCnt=0
    ILambda_29_0Cnt=0
    ILambda_29_0FinalResCnt=0
    ILambda_29_0ItersCnt=0
    ILambda_3_0Cnt=0
    ILambda_30_0Cnt=0
    ILambda_30_0FinalResCnt=0
    ILambda_30_0ItersCnt=0
    ILambda_3_0FinalResCnt=0
    ILambda_3_0ItersCnt=0
    ILambda_31_0Cnt=0
    ILambda_31_0FinalResCnt=0
    ILambda_31_0ItersCnt=0
    ILambda_32_0Cnt=0
    ILambda_32_0FinalResCnt=0
    ILambda_32_0ItersCnt=0
    ILambda_33_0Cnt=0
    ILambda_33_0FinalResCnt=0
    ILambda_33_0ItersCnt=0
    ILambda_34_0Cnt=0
    ILambda_34_0FinalResCnt=0
    ILambda_34_0ItersCnt=0
    ILambda_35_0Cnt=0
    ILambda_35_0FinalResCnt=0
    ILambda_35_0ItersCnt=0
    ILambda_36_0Cnt=0
    ILambda_36_0FinalResCnt=0
    ILambda_36_0ItersCnt=0
    ILambda_37_0Cnt=0
    ILambda_37_0FinalResCnt=0
    ILambda_37_0ItersCnt=0
    ILambda_38_0Cnt=0
    ILambda_38_0FinalResCnt=0
    ILambda_38_0ItersCnt=0
    ILambda_39_0Cnt=0
    ILambda_39_0FinalResCnt=0
    ILambda_39_0ItersCnt=0
    ILambda_4_0Cnt=0
    ILambda_40_0Cnt=0
    ILambda_40_0FinalResCnt=0
    ILambda_40_0ItersCnt=0
    ILambda_4_0FinalResCnt=0
    ILambda_4_0ItersCnt=0
    ILambda_41_0Cnt=0
    ILambda_41_0FinalResCnt=0
    ILambda_41_0ItersCnt=0
    ILambda_42_0Cnt=0
    ILambda_42_0FinalResCnt=0
    ILambda_42_0ItersCnt=0
    ILambda_43_0Cnt=0
    ILambda_43_0FinalResCnt=0
    ILambda_43_0ItersCnt=0
    ILambda_44_0Cnt=0
    ILambda_44_0FinalResCnt=0
    ILambda_44_0ItersCnt=0
    ILambda_45_0Cnt=0
    ILambda_45_0FinalResCnt=0
    ILambda_45_0ItersCnt=0
    ILambda_46_0Cnt=0
    ILambda_46_0FinalResCnt=0
    ILambda_46_0ItersCnt=0
    ILambda_47_0Cnt=0
    ILambda_47_0FinalResCnt=0
    ILambda_47_0ItersCnt=0
    ILambda_48_0Cnt=0
    ILambda_48_0FinalResCnt=0
    ILambda_48_0ItersCnt=0
    ILambda_49_0Cnt=0
    ILambda_49_0FinalResCnt=0
    ILambda_49_0ItersCnt=0
    ILambda_5_0Cnt=0
    ILambda_50_0Cnt=0
    ILambda_50_0FinalResCnt=0
    ILambda_50_0ItersCnt=0
    ILambda_5_0FinalResCnt=0
    ILambda_5_0ItersCnt=0
    ILambda_51_0Cnt=0
    ILambda_51_0FinalResCnt=0
    ILambda_51_0ItersCnt=0
    ILambda_52_0Cnt=0
    ILambda_52_0FinalResCnt=0
    ILambda_52_0ItersCnt=0
    ILambda_53_0Cnt=0
    ILambda_53_0FinalResCnt=0
    ILambda_53_0ItersCnt=0
    ILambda_54_0Cnt=0
    ILambda_54_0FinalResCnt=0
    ILambda_54_0ItersCnt=0
    ILambda_55_0Cnt=0
    ILambda_55_0FinalResCnt=0
    ILambda_55_0ItersCnt=0
    ILambda_56_0Cnt=0
    ILambda_56_0FinalResCnt=0
    ILambda_56_0ItersCnt=0
    ILambda_57_0Cnt=0
    ILambda_57_0FinalResCnt=0
    ILambda_57_0ItersCnt=0
    ILambda_58_0Cnt=0
    ILambda_58_0FinalResCnt=0
    ILambda_58_0ItersCnt=0
    ILambda_59_0Cnt=0
    ILambda_59_0FinalResCnt=0
    ILambda_59_0ItersCnt=0
    ILambda_6_0Cnt=0
    ILambda_60_0Cnt=0
    ILambda_60_0FinalResCnt=0
    ILambda_60_0ItersCnt=0
    ILambda_6_0FinalResCnt=0
    ILambda_6_0ItersCnt=0
    ILambda_61_0Cnt=0
    ILambda_61_0FinalResCnt=0
    ILambda_61_0ItersCnt=0
    ILambda_62_0Cnt=0
    ILambda_62_0FinalResCnt=0
    ILambda_62_0ItersCnt=0
    ILambda_63_0Cnt=0
    ILambda_63_0FinalResCnt=0
    ILambda_63_0ItersCnt=0
    ILambda_7_0Cnt=0
    ILambda_7_0FinalResCnt=0
    ILambda_7_0ItersCnt=0
    ILambda_8_0Cnt=0
    ILambda_8_0FinalResCnt=0
    ILambda_8_0ItersCnt=0
    ILambda_9_0Cnt=0
    ILambda_9_0FinalResCnt=0
    ILambda_9_0ItersCnt=0
    SeparatorCnt=0
    TimeCnt=0
    # Reset counters for 'Solving for ...'
    for (varName in subIter)
    {
        subIter[varName]=0
    }
}

# Extract value after columnSel
function extract(inLine,columnSel,outVar,a,b)
{
    a=index(inLine, columnSel)
    b=length(columnSel)
    split(substr(inLine, a+b),outVar)
    gsub("[,:]","",outVar[1])
}

# Iteration separator (increments 'Iteration')
/^[ \t]*Time = / {
    Iteration++
    resetCounters()
}

# Time extraction (sets 'Time')
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    Time=val[1]
}

# Skip whole line with singularity variable
/solution singularity/ {
    next;
}

# Extract: 'Solving for ...'
/Solving for/ {
    extract($0, "Solving for ", varNameVal)

    varName=varNameVal[1]
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "Initial residual = ", val)
    print Time "\t" val[1] > file

    varName=varNameVal[1] "FinalRes"
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "Final residual = ", val)
    print Time "\t" val[1] > file

    varName=varNameVal[1] "Iters"
    file=varName "_" subIter[varName]++
    file="logs/" file
    extract($0, "No Iterations ", val)
    print Time "\t" val[1] > file
}

# Extract: 'clockTime'
/ClockTime = / {
    extract($0, "ClockTime =", val)
    file="logs/clockTime_" clockTimeCnt
    print Time "\t" val[1] > file
    clockTimeCnt++
}

# Extract: 'executionTime'
/ExecutionTime = / {
    extract($0, "ExecutionTime = ", val)
    file="logs/executionTime_" executionTimeCnt
    print Time "\t" val[1] > file
    executionTimeCnt++
}

# Extract: 'Separator'
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    file="logs/Separator_" SeparatorCnt
    print Time "\t" val[1] > file
    SeparatorCnt++
}

# Extract: 'Time'
/^[ \t]*Time = / {
    extract($0, "Time = ", val)
    file="logs/Time_" TimeCnt
    print Time "\t" val[1] > file
    TimeCnt++
}

# End
