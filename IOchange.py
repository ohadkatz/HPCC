import sys
import fileinput
def ChangeInput():
    with open("/home/ohadkatz/HPCC/Results.txt", "r+") as OutputFile:
        NewOut= OutputFile.readlines()
        Values = sys.argv[1:]
        Processors = Values[0]
        Threads = Values[1]
        Type = Values[2]
        print(Type)
        Start = int(Values[3])
        End = int(Values[4])
        NewOut[0]= "Type, Proc, Threads, "+ NewOut[0]
        for i in range(Start+1,End):
            if "Done" in NewOut[i]:
                continue
            NewOut[i]= Type + "," + Processors + "," + Threads + "," + NewOut[i]
    with open("/home/ohadkatz/HPCC/FinalResults.txt", "a+") as f:
        for i in NewOut:
            f.writelines(i)

        
if __name__=="__main__":
    ChangeInput()