import numpy



def main():

    weird_hera_lst = 64829.1234567
    fractional_date = weird_hera_lst % 1
    print(weird_hera_lst)
    print(fractional_date)
    hourangle = fractional_date * 24
    print("Hour Angle: ")
    print(hourangle)

if __name__=="__main__":
    main()
