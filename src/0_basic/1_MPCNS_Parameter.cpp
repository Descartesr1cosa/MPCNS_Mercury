#include "0_basic/1_MPCNS_Parameter.h"
#include "0_basic/CONST_DEFINE.h"

void Param::ReadParam(int32_t _myid)
{
    AddParam("myid", _myid);
    if (GetInt("myid") == 0)
    {
        std::cout << "====================================="
                  << "===============================================================\n";
        std::cout << "\t \t Welcome to the Multi-Physics Coupling Numerical Simulation program\n";
        std::cout << "====================================="
                  << "===============================================================\n\n\n";
        std::cout << "---->Reading the Control Parameters...\n";
    }

    /*创建一个map用于判断输入的参数是哪个类型*/
    std::map<std::string, int> keywordmap;
    keywordmap.insert(std::pair<std::string, int>("string", STRING));
    keywordmap.insert(std::pair<std::string, int>("int", INT));
    keywordmap.insert(std::pair<std::string, int>("double", DOUBLE));
    keywordmap.insert(std::pair<std::string, int>("bool", BOOL));
    keywordmap.insert(std::pair<std::string, int>("IF", IF));
    keywordmap.insert(std::pair<std::string, int>("ELSEIF", ELSEIF));
    keywordmap.insert(std::pair<std::string, int>("ELSE", ELSE));
    keywordmap.insert(std::pair<std::string, int>("ENDIF", ENDIF));
    keywordmap.insert(std::pair<std::string, int>("END_OF_FILE", END_OF_FILE));
    keywordmap.insert(std::pair<std::string, int>("LIST", LIST));

    /*sep(seperator)用于分隔一个参数 类型/名称/值*/
    std::string sep = "=\r\n\t#$,;\" ()";

    //---------------------------------------------------------------------------------------------
    /*read filenames to know which file has to be read*/
    std::fstream file;
    file.open("./CASE/setup/filenames", std::ios_base::in);
    if (!file)
    {
        std::cout << "can't open file : \t"
                  << "./CASE/setup/filenames" << std::endl;
        return;
    }
    while (!file.eof())
    {
        std::string filename;
        file >> filename; // getline(file, filename);
        if (filename == "###END###")
            break;
        filename = "./CASE/setup/" + filename;
        input_file.push_back(filename);
    }
    file.close();
    /*read filenames to know which file has to be read*/
    //---------------------------------------------------------------------------------------------

    /*read input_file of General_Setting*/
    for (int i = 0; i < input_file.size(); i++)
        ReadFile(input_file[i], keywordmap, sep);
    // 参数预处理
    Pre_Process();
    // 参数读入完成！
    if (GetInt("myid") == 0)
        std::cout << "****************Finish the General reading Process! !*****************\n"
                  << std::endl;
}

std::string Param::FindNextWord(std::string &source, std::string &word, std::string &sep)
{
    std::string::size_type firstindex, nextindex, npos = std::string::npos;
    std::string nullstring = "";
    firstindex = source.find_first_not_of(sep, 0);

    if (firstindex == npos)
    {
        word = nullstring;
        return nullstring;
    }

    nextindex = source.find_first_of(sep, firstindex);
    if (nextindex == npos)
    {
        word = source.substr(firstindex);
        return nullstring;
    }
    else
    {
        word = source.substr(firstindex, nextindex - firstindex);
        return source.substr(nextindex);
    }
}

bool Param::ReadIFModule_Justify(std::string &filename, std::string &line, std::map<std::string, int> &keywordmap, std::string &sep)
{
    std::string keyword, name, word, value;
    // 提取用于判断条件的量是什么类型，这里只支持int bool string类型
    line = FindNextWord(line, keyword, sep);
    if (keywordmap.find(keyword) == keywordmap.end())
    {
        std::cout << "error in inputfile of IF Module" << filename << "\t" << keyword << std::endl;
        exit(-1);
    }
    int condition_i, input_i;
    std::string condition_s, input_s;
    bool condition_b, input_b, condition;
    input_b = false;
    condition = false;

    if (keywordmap[keyword] == INT)
    {
        line = FindNextWord(line, word, sep);
        condition_i = GetInt(word);
        line = FindNextWord(line, value, sep);
        input_i = stoi(value);
        if (condition_i == input_i)
            condition = true;
    }
    else if (keywordmap[keyword] == BOOL)
    {
        line = FindNextWord(line, word, sep);
        condition_b = GetBoo(word);
        line = FindNextWord(line, value, sep);
        if (value != "0")
            input_b = true;
        if (condition_b == input_b)
            condition = true;
    }
    else if (keywordmap[keyword] == STRING)
    {
        line = FindNextWord(line, word, sep);
        condition_s = GetStr(word);
        line = FindNextWord(line, value, sep);
        input_s = value;
        if (condition_s == input_s)
            condition = true;
    }
    return condition;
}

void Param::ReadIFModule(std::string &filename, std::fstream &file, std::string &line, std::map<std::string, int> &keywordmap, std::string &sep)
{
    std::string keyword, name, word, value;
    // 提取用于判断条件的量是什么类型，这里只支持int bool string类型
    bool condition;
    condition = ReadIFModule_Justify(filename, line, keywordmap, sep);
    for (; !condition;)
    {
        getline(file, line);

        // 空行直接过
        if (line.find_first_not_of(sep) == std::string::npos)
            continue;

        // 这一行开始有“//”或“#”为注释行，不读，直接跳下一行
        if ((line.substr(0, 2) == "//") || (line.substr(0, 1) == "#"))
            continue;

        line = FindNextWord(line, keyword, sep);
        if (keywordmap.find(keyword) == keywordmap.end())
        {
            std::cout << "error in inputfile of IF Module 2" << filename << std::endl;
        }

        if (keywordmap[keyword] == ELSEIF) // ELSEIF
        {
            condition = ReadIFModule_Justify(filename, line, keywordmap, sep);
        }
        else if (keywordmap[keyword] == ELSE) // ELSE
        {
            break;
        }
        else if (keywordmap[keyword] == ENDIF)
        {
            return;
        }
    }
    // 若IF判断正确，直接运行到此，不用上面的循环
    // 若IF判断不正确，循环找ELSE ELSEIF ENDIF
    //       ELSE 在循环中说明IF判断以及其他的ELSEIF均不正确，直接读如ELSE后面的内容到ENDIF
    //       ELSEIF 判断更新条件即可，若符合条件就退出到此读入到下一个ELSE ELSEIF  ENDIF内容，否则继续循环
    //       ENDIF 直接return,不会到此
    // 此处应当写一直读入到下一个ELSE ELSEIF ENDIF的内容
    // 为了后面完整，应当一直getline读到ENDIF，但是不AddParam写入数据
    condition = true;
    bool end_read_condition = false; // 是否结束AddParam写入读入数据
    do
    {
        getline(file, line);

        // 空行直接过
        if (line.find_first_not_of(sep) == std::string::npos)
            continue;

        // 这一行开始有“//”或“#”为注释行，不读，直接跳下一行
        if ((line.substr(0, 2) == "//") || (line.substr(0, 1) == "#"))
            continue;

        line = FindNextWord(line, keyword, sep);
        if (keywordmap.find(keyword) == keywordmap.end())
        {
            std::cout << "error in inputfile of IF Module 3" << filename << std::endl;
            exit(-1);
        }

        line = FindNextWord(line, name, sep);
        line = FindNextWord(line, value, sep);

        // 若end_read_condition=false表示还需要AddParam写入数据
        // 直到碰到ELSE ELSEIF ENDIF才停止 end_read_condition = true
        if (!end_read_condition)
        {
            switch (keywordmap[keyword])
            {
            case INT:
            {
                int temp_i = stoi(value);
                AddParam(name, temp_i);
            }
            break;

            case DOUBLE:
            {
                double temp_d = stod(value);
                AddParam(name, temp_d);
            }
            break;

            case BOOL:
            {
                bool temp_b = false;
                if (value != "0")
                    temp_b = true;
                AddParam(name, temp_b);
            }
            break;

            case STRING:
                AddParam(name, value);
                break;

            case LIST:
                ReadLISTModule(filename, file, line, keywordmap, sep);
                break;

            case ENDIF:
                // ENDIF表示已经结束，可以直接结束循环
                condition = false;
                end_read_condition = true;
                break;

            default:
                end_read_condition = true;
                break;
            }
        }
        else // 若end_read_condition=true表示不需要AddParam写入数据，直到ENDIF
        {
            switch (keywordmap[keyword])
            {
            case ENDIF:
                // ENDIF表示已经结束，可以直接结束循环
                condition = false;
                break;
            default:
                break;
            }
        }

    } while (condition);
}

void Param::ReadLISTModule(std::string &filename, std::fstream &file, std::string &line, std::map<std::string, int> &keywordmap, std::string &sep)
{
    std::string keyword, name, word, value, line_name, List_name;

    getline(file, line); // 本行说明List的名称和List的个数
    line = FindNextWord(line, List_name, sep);
    line = FindNextWord(line, value, sep);

    getline(file, line); // 本行说明类型，仅仅支持int double类型
    line = FindNextWord(line, keyword, sep);
    if (keywordmap[keyword] == INT)
    {
        List<int32_t> temp(stoi(value));
        getline(file, line); // 本行说明变量
        line_name = line;
        getline(file, line); // 本行为变量的值

        for (int i = 0; i < temp.num; i++)
        {
            line_name = FindNextWord(line_name, name, sep);
            line = FindNextWord(line, value, sep);
            temp.data[name] = stoi(value);
        }
        getline(file, line); // 本行读入List
        AddParam(List_name, temp);
        return;
    }
    else if (keywordmap[keyword] == DOUBLE)
    {
        List<double> temp(stoi(value));
        getline(file, line); // 本行说明变量
        line_name = line;
        getline(file, line); // 本行为变量的值

        for (int i = 0; i < temp.num; i++)
        {
            line_name = FindNextWord(line_name, name, sep);
            line = FindNextWord(line, value, sep);
            temp.data[name] = stod(value);
        }
        getline(file, line); // 本行读入List
        AddParam(List_name, temp);
        return;
    }
}

void Param::ReadFile(std::string &filename, std::map<std::string, int> keywordmap, std::string sep)
{
    std::string line, keyword, name, word, value;
    std::fstream file;
    bool endfile = false;
    file.open(filename.c_str(), std::ios_base::in);
    if (!file)
    {
        std::cout << "can't open file : " << filename << "_" << std::endl;
        return;
    }
    while (!file.eof())
    {
        getline(file, line);
        // 空行直接过
        if (line.find_first_not_of(sep) == std::string::npos)
            continue;

        // 这一行开始有“//”或“#”为注释行，不读，直接跳下一行
        if ((line.substr(0, 2) == "//") || (line.substr(0, 1) == "#"))
            continue;

        /*读入第一个参数，表征类型除了int double bool string 还有IF ELSE ELSEIF ENDIF END_OF_FILE*/
        // Remark: 由于IF模块化读入，因此不应当出现ELSE ELSEIF ENDIF
        line = FindNextWord(line, keyword, sep);
        if (keywordmap.find(keyword) == keywordmap.end())
        {
            std::cout << "error in inputfile of ReadFile 1" << filename << std::endl;
            exit(-1);
        }

        switch (keywordmap[keyword])
        {
        case INT:
        {
            /*读输入参数名称和参数值*/
            line = FindNextWord(line, name, sep);
            line = FindNextWord(line, value, sep);
            int temp = stoi(value);
            AddParam(name, temp);
        }
        break;

        case DOUBLE:
        {
            /*读输入参数名称和参数值*/
            line = FindNextWord(line, name, sep);
            line = FindNextWord(line, value, sep);
            double temp = stod(value);
            AddParam(name, temp);
        }
        break;

        case BOOL:
        {
            /*读输入参数名称和参数值*/
            line = FindNextWord(line, name, sep);
            line = FindNextWord(line, value, sep);
            bool temp = false;
            if (value != "0")
                temp = true;
            AddParam(name, temp);
        }
        break;

        case LIST:
            /* LIST */
            ReadLISTModule(filename, file, line, keywordmap, sep);
            break;

        case IF:
            /* IF */
            ReadIFModule(filename, file, line, keywordmap, sep);
            break;

        case STRING:
            /*读输入参数名称和参数值*/
            line = FindNextWord(line, name, sep);
            line = FindNextWord(line, value, sep);
            AddParam(name, value);
            break;

        case END_OF_FILE:
            endfile = true;
            break;

        default:
            std::cout << "error in inputfile of ReadFile2" << filename << std::endl;
            exit(-1);
            break;
        }
        if (endfile)
            break;
    }
    file.close();
    if (this->GetInt("myid") == 0)
        std::cout << "\t-->Finish reading :\t" << filename << std::endl;
}

void Param::Pre_Process()
{
    // if (IntParam["ngg"] <= 2)
    // {
    //     AddParam("ngg", 2);
    // }
    // else if (IntParam["ngg"] <= 3)
    // {
    //     AddParam("ngg", 3);
    // }
    // else if (IntParam["ngg"] <= 5)
    // {
    //     AddParam("ngg", 5);
    // }
    // else
    // {
    //     AddParam("ngg", 7);
    // }
}
