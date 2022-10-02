#include "ooqp/QpGenData.h"
#include "ooqp/QpGenVars.h"
#include "ooqp/QpGenResiduals.h"
#include "ooqp/GondzioSolver.h"
#include "ooqp/QpGenSparseMa27.h"

#include <iostream>
#include <string.h>

using namespace std;

const int nx = 2;//变量的长度
double c[] = {1.5,-2};//目标函数中的线性部分，长度为nx

//使用OOQP求解问题的步骤总结：
//1.调用类初始化一个问题求解器，不同的类决定了不同的数据存储方式和不同的线性系统的求解方式。
//2.利用对象调用其中的 makeData 方法来创建向量和矩阵来存储和问题有关的数据。
//3.利用对象调用其中的 makeVariables 方法来存储问题中的变量。
//4.利用对象调用其中的 makeResiduals 方法来存储在不同给定点下不同的优化条件下的残差。
//5.创建求解器来求解 QP 问题 GondzioSolver *s = new GondzioSolver(qp,prob);
//6.调用上面一步创建的 s 来调用 monitorSelf 方法把解算得到的信息打印到屏幕上。
//7.调用 s->solve(prob,vars,resid) 来得到结果，其返回值为0时说明能够计算得到值，如果返回值非零则问题不可解，返回值为负数说明有错误出现例如内存错误等等。如果返回值为正数会存储在 ierr 中。

//稀疏矩阵由三个数据结构来表示-两个整数向量和一个双精度向量，他们的长度都相同。对于矩阵C，这些数据结构为 irowC,jcolC和dC，稀疏矩阵中的非零元素总数为 nnzC。C的第k个非零元素出现在 irowC[k] 和 icolC[k] 处，其值为 dC[k]。注意：行和列时从零开始编号的。
//对称矩阵也用三个向量来表示，在 irowQ,jcolQ和dQ中指定矩阵下三角（左下角）的元素。
//上下界采用两个向量来表示，例如，如果x的下界的第k个元素有下界，ixlow[k]设置为1，xlow[k]为其下界值；否则都设置为0。
double xupp[] = {20,0};
char ixupp[] = {1,0};

double xlow[] = {0,0};
char ixlow[] = {1,1};

const int nnzQ = 3;//对于矩阵Q，仅仅标明左下三角的元素即可，包含零元素
int irowQ[] = {0,1,1};
int jcolQ[] = {0,0,1};
double dQ[] = {8,2,10};

int my = 0;//线性等式约束的个数
double *b = 0;//等式约束右侧的值
int nnzA = 0;//矩阵A中非零元素的个数
int *irowA = 0;//非零元素所在的行数，行数从零开始
int *jcolA = 0;//非零元素所在的列数，列数从零开始
double *dA = 0;//非零元素的值

// int my = 1;//线性等式约束的个数
// double b[] = {0};
// int nnzA = 2;
// int irowA[] = {0,0};
// int jcolA[] = {0,1};
// double dA[] = {1,1};

// const int mz = 0;
// double *clow = 0;
// char *iclow = 0;
// double *cupp = 0;
// char *icupp = 0;

// const int nnzC = 0;
// int *irowC = 0;
// int *jcolC = 0;
// double *dC = 0;

const int mz = 2;
double clow[] = {2,0};
char iclow[] = {1,0};
double cupp[] = {0,6};
char icupp[] = {0,1};

const int nnzC = 4;//矩阵C中的非零元素的个数
int irowC[] = {0,0,1,1};//非零元素所在行
int jcolC[] = {0,1,0,1};//非零元素所在列
double dC[] = {2,1,-1,2};//非零元素的值

int main(int argc,char * argv[])
{
	int usage_ok = 1;
	int quiet = 0;
	if(argc > 2)
	{
		usage_ok = 0;
	}
	if(argc == 2)
	{
		if(0 == strcmp("--quiet",argv[1]))
		{
			quiet = 1;
		}
		else
		{
			usage_ok = 0;
		}
	}	

	if(!usage_ok)
	{
		cerr<<"Usage:"<<argv[0]<<"[--quiet]"<<endl;
		return 1;
	}
	QpGenSparseMa27 * qp = new QpGenSparseMa27(nx,my,mz,nnzQ,nnzA,nnzC);
	QpGenData * prob = (QpGenData *)qp -> copyDataFromSparseTriple(
			c,irowQ,nnzQ,jcolQ,dQ,
			xlow,ixlow,xupp,ixupp,
			irowA,nnzA,jcolA,dA,b,
			irowC,nnzC,jcolC,dC,
			clow,iclow,cupp,icupp);
	QpGenVars * vars = (QpGenVars *)qp -> makeVariables(prob);
	QpGenResiduals *resid = (QpGenResiduals *)qp -> makeResiduals(prob);

	GondzioSolver *s = new GondzioSolver(qp,prob);
	
	if(!quiet)
	{
		cout<<"quiet = "<<quiet<<endl;
		s->monitorSelf();
	}
	int ierr = s->solve(prob,vars,resid);

	if(ierr == 0)
	{
		cout.precision(4);
		cout<<"Solution:"<<endl;
		vars->x->writefToStream(cout,"x[%{index}] = %{value}");
	}
	else
	{
		cout<<"Could not solve the problem!"<<endl;
	}

	return ierr;
}
