/*****************************************************************************
@Copyright: NLCFD
@File name: 1_MPCNS_Parameter.h
@Author: WangYH Descartes
@Version: 1.0
@Date: 2022年09月11日
@Description:	基于自定义的可变长度Array数组定义了控制参数并形成了读入参数的方法
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了文件名称，基本的功能头文件名用数字1+全部首字母大写的方式
*****************************************************************************/
#pragma once
#include <map>
#include <fstream>
#include "0_basic/ARRAY_DEFINE.h"

/**
 * @brief			存储以LIST方式读入的计算参数
 *					求解器设置、参考量、计算用到的常数等
 *
 * @param num		参数的个数
 * @param data	    map类型，string和int/double组成，存储输入的参数
 */
template <typename T>
class List
{
public:
    int32_t num;
    std::map<std::string, T> data;

public:
    List() {};
    List(int32_t N)
    {
        num = N;
    };
    ~List() = default;
};

/**
 * @brief			存储计算过程中的所有计算参数
 *					求解器设置、参考量、计算用到的常数等
 *
 * @param input_file		部分输入文件名 .\setup\***.txt
 * @param StrParam	        map类型，string和string组成，存储输入的字符串变量
 * @param IntParam	        map类型，string和int组成，存储输入的整型变量
 * @param DouParam	        map类型，string和double组成，存储输入的双精度浮点型变量
 * @param BooParam	        map类型，string和bool组成，存储输入的布尔类型变量
 */
class Param
{
public:
    Param() {};
    ~Param() = default;

private:
    // 输入文件名 .\setup\*** .txt
    std::vector<std::string> input_file;

private:
    // 存储string类型参数的map
    std::map<std::string, std::string> StrParam;
    // 存储int类型参数的map
    std::map<std::string, int> IntParam;
    // 存储double类型参数的map
    std::map<std::string, double> DouParam;
    // 存储bool类型参数的map
    std::map<std::string, bool> BooParam;
    // 存储int List类型参数的map
    std::map<std::string, List<int>> List_IntParam;
    // 存储double List类型参数的map
    std::map<std::string, List<double>> List_DouParam;

public:
    /**
     * @brief	更新参数，重载4次，分别对应四种参数
     * @param name_in 通过字符串索引参数，再对参数更新
     */
    void AddParam(std::string name_in, std::string data_in) { StrParam[name_in] = data_in; };       // 更新string类型参数
    void AddParam(std::string name_in, int data_in) { IntParam[name_in] = data_in; };               // 更新int类型参数
    void AddParam(std::string name_in, double data_in) { DouParam[name_in] = data_in; };            // 更新double类型参数
    void AddParam(std::string name_in, bool data_in) { BooParam[name_in] = data_in; };              // 更新bool类型参数
    void AddParam(std::string name_in, List<int> data_in) { List_IntParam[name_in] = data_in; };    // 更新int List类型参数
    void AddParam(std::string name_in, List<double> data_in) { List_DouParam[name_in] = data_in; }; // 更新double List类型参数

    /**
     * @brief	索引参数，有四个函数，分别对应四种参数
     * @param name_in 通过字符串索引参数
     * @return	返回索引的参数值
     */
    std::string GetStr(std::string name_in)
    {
        if (StrParam.count(name_in) == 0)
        {
            std::cout << "#Fatal Error: No such String as:\t" << name_in << std::endl;
            exit(-1);
        }
        return StrParam[name_in];
    }; // 索引string类型参数
    int GetInt(std::string name_in)
    {
        if (IntParam.count(name_in) == 0)
        {
            std::cout << "#Fatal Error: No such Int as:\t" << name_in << std::endl;
            exit(-1);
        }
        return IntParam[name_in];
    }; // 索引Int类型参数
    double GetDou(std::string name_in)
    {
        if (DouParam.count(name_in) == 0)
        {
            std::cout << "#Fatal Error: No such Double as:\t" << name_in << std::endl;
            exit(-1);
        }
        return DouParam[name_in];
    }; // 索引double类型参数
    bool GetBoo(std::string name_in)
    {
        if (BooParam.count(name_in) == 0)
        {
            std::cout << "#Fatal Error: No such Bool as:\t" << name_in << std::endl;
            exit(-1);
        }
        return BooParam[name_in];
    }; // 索引bool类型参数
    List<int> GetInt_List(std::string name_in)
    {
        if (List_IntParam.count(name_in) == 0)
        {
            std::cout << "#Fatal Error: No such List(Int) as:\t" << name_in << std::endl;
            exit(-1);
        }
        return List_IntParam[name_in];
    }; // 索引int List类型参数
    List<double> GetDou_List(std::string name_in)
    {
        if (List_DouParam.count(name_in) == 0)
        {
            std::cout << "#Fatal Error: No such List(Double) as:\t" << name_in << std::endl;
            exit(-1);
        }
        return List_DouParam[name_in];
    }; // 索引double List类型参数

    /**
     * @brief	读输入文件
     *          分为单个参数读入，IF ENDIF模块读入 LIST模块读入三类
     *          IF ENDIF中可以有单个参数读入和LIST模块读入
     */
    void ReadParam(int32_t _myid);

private:
    /**
     * @brief	通过关键字读取文件中的输入参数
     */
    void ReadFile(std::string &filename, std::map<std::string, int> keywordmap, std::string sep);

    //---------------------------------------------------------------------------------------------
    /**
     * @brief	通过关键字读取文件中的IF模块
     */
    void ReadIFModule(std::string &filename, std::fstream &file, std::string &line, std::map<std::string, int> &keywordmap, std::string &sep);

    /**
     * @brief	通过本行给出读取行的判断结果
     */
    bool ReadIFModule_Justify(std::string &filename, std::string &line, std::map<std::string, int> &keywordmap, std::string &sep);
    //---------------------------------------------------------------------------------------------

    /**
     * @brief	通过关键字读取文件中的LIST模块
     */
    void ReadLISTModule(std::string &filename, std::fstream &file, std::string &line, std::map<std::string, int> &keywordmap, std::string &sep);

    /**
     * @brief	从source字符串中找到第一个不出现在sep字符串中的字符，赋给word
     */
    std::string FindNextWord(std::string &source, std::string &word, std::string &sep);

    void Pre_Process();
};
