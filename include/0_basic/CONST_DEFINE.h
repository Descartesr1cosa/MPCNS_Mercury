/*****************************************************************************
@Copyright: NLCFD
@File name: 0_CONST_DEFINE.h
@Author: WangYH Descartes
@Version: 1.0
@Date: 2022年10月01日
@Description:	定义常数、宏变量。
                注意：良好的编程习惯，宏变量名称最好大写
@History:		（修改历史记录列表， 每条修改记录应包括修改日期、修改者及修改内容简述。）
                1、修改了文件名称，最基本的功能头文件名用数字0+全部大写的方式
*****************************************************************************/

#pragma once

// 枚举，给各种类型编号
constexpr int STRING = 1;      // 字符串类型编号
constexpr int INT = 2;         // 整形类型编号
constexpr int DOUBLE = 3;      // 双精度浮点数类型编号
constexpr int BOOL = 4;        // 布尔类型编号
constexpr int IF = 5;          // 判断型
constexpr int ELSEIF = 6;      // 判断型
constexpr int ELSE = 7;        // 判断型
constexpr int ENDIF = 8;       // 判断型
constexpr int END_OF_FILE = 9; // 文件结束
constexpr int LIST = 10;       // 文件结束