#include "/home/zjs/minimum_snap_ws/src/optimize_trim/include/optimize_trim/OoqpOptimize.h"

OptimizeOoqp::OptimizeOoqp()
{
    EqualConstrain = false;
    UnequalConstrain = false;
}

void OptimizeOoqp::SingleSolveCurve(const double T,const int EqualConstrainNum,const int UnequalConNum,
                                                                    const double ps,const double vs,const double as,
                                                                    const double pe,const double ve,const double ae)
{
    const int nx = 6;
    double c[] = {0,0,0,0,0,0};

    double xupp[] = {0,0,0,0,0,0};
    char ixupp[] = {0,0,0,0,0,0};
    double xlow[] = {0,0,0,0,0,0};
    char ixlow[] = {0,0,0,0,0,0};

    const int nnzQ = 21;
    int irowQ[] = {0,1,2,3,4,5,1,2,3,4,5,2,3,4,5,3,4,5,4,5,5};
    int jcolQ[] = {0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,3,3,3,4,4,5};
    double dQ[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,576 * std::pow(T,1),1440 * std::pow(T,2),4800 * std::pow(T,3)};

    int my = EqualConstrainNum;
    double b[] = {ps,vs,as,pe,ve,ae};
    int nnzA = 18;
    int irowA[] = {0,1,2,3,3,3,3,3,3,4,4,4,4,4,5,5,5,5};
    int jcolA[] = {5,4,3,0,1,2,3,4,5,0,1,2,3,4,0,1,2,3};
    double dA[] = {1,1,2,pow(T,5),pow(T,4),pow(T,3),pow(T,2),pow(T,1),1,5*pow(T,4),4*pow(T,3),3*pow(T,2),2*pow(T,1),1,20*pow(T,3),12*pow(T,2),6*pow(T,1),2};

    int mz = UnequalConNum;
    double *clow = 0;
    char *iclow = 0;
    double *cupp = 0;
    char *icupp = 0;

    const int nnzC = 0;
    int *irowC = 0;
    int *jcolC = 0;
    double *dC = 0;

    QpGenSparseMa27 * qp = new QpGenSparseMa27(nx,my,mz,nnzQ,nnzA,nnzC);
    QpGenData *prob = (QpGenData *)qp ->copyDataFromSparseTriple(
                                                                                                                                                c,irowQ,nnzQ,jcolQ,dQ,
                                                                                                                                                xlow,ixlow,xupp,ixupp,
                                                                                                                                                irowA,nnzA,jcolA,dA,b,
                                                                                                                                                irowC,nnzC,jcolC,dC,
                                                                                                                                                clow,iclow,cupp,icupp);
    QpGenVars * vars = (QpGenVars *)qp ->makeVariables(prob);
    QpGenResiduals *resid = (QpGenResiduals *)qp -> makeResiduals(prob);
    GondzioSolver *s = new GondzioSolver(qp,prob);

    s->monitorSelf();
    int ierr = s -> solve(prob,vars,resid);
    if(ierr == 0)
    {
        cout.precision(4);
        cout<<"Solution: "<<endl;
        vars->x->writefToStream(cout,"x[%{index}] = %{value}");
    }
    else
    {
        cout<<"Could not solve the problem !"<<endl;
    }
}

void OptimizeOoqp::Test_Solve()
{
    const int nx = 2;//变量的长度
    double c[] = {1.5,-2};//目标函数中的线性部分，长度为nx

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

    const int mz = 2;
    double clow[] = {2,0};
    char iclow[] = {1,0};
    double cupp[] = {0,6};
    char icupp[] = {0,1};

    const int nnzC = 4;//矩阵C中的非零元素的个数
    int irowC[] = {0,0,1,1};//非零元素所在行
    int jcolC[] = {0,1,0,1};//非零元素所在列
    double dC[] = {2,1,-1,2};//非零元素的值

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
    
    s->monitorSelf();
    int ierr = s->solve(prob,vars,resid);

	if(ierr == 0)
	{
		cout.precision(4);
		cout<<"Solution:"<<endl;
		vars->x->writefToStream(cout,"x[%{index}] = %{value}");
	}
}

void OptimizeOoqp::SolveCurve(const double sumT,const int EquConNum,const int UnequConNum,const Eigen::Vector3d& begin,const Eigen::Vector3d& end,const Eigen::VectorXd& midle)
{
    Eigen::VectorXd timePor;
    int num = midle.size() + 2 - 1;
    int sum = 0;

    timePor.resize(num);
    for(int i = 0;i < midle.size();i++)
    {
        if(i == 0)
        {
            sum = sum + std::abs((midle(0) - begin(0)));
        }
        if(i > 0 && i < midle.size())
        {
            sum = sum + std::abs((midle(i) - midle(i-1)));
        }
    }
    sum = sum + std::abs(end(0) - midle(midle.size() - 1));

    for (int i = 0; i < midle.size(); i++)
    {
        if(i == 0)
        {
            timePor(i) = sumT * (std::abs((midle(0) - begin(0))) / sum);
        }
        if(i > 0 && i < midle.size())
        {
            timePor(i) = sumT * (std::abs((midle(i) - midle(i-1))) / sum);
        }
    }
    timePor(num-1) = sumT * (std::abs(end(0) - midle(midle.size() - 1)) / sum);

    std::cout<<"sum = "<<sum<<endl;
    for (int i = 0; i < timePor.size(); i++)
    {
        cout<<"timePor["<<i<<"] = "<<timePor(i)<<endl;
    }   
}

//maxorder is 5,so n = 5.every curve includes 6 parameters(p5,p4,p3,p2,p1,p0).
void OptimizeOoqp::ConEqualSolve(const Eigen::VectorXd& T,const int curveNums,const Eigen::Vector3d& begin,const Eigen::Vector3d& end,const Eigen::VectorXd& midlep)
{
    const int nx = 6 * curveNums;
    double c[nx] = {0};

    double xupp[nx] = {0};
    char ixupp[nx] = {0};
    double xlow[nx] = {0};
    char ixlow[nx] = {0};

    const int nnzQ = 21 * curveNums;
    int irowQ[nnzQ] = {0};
    int jcolQ[nnzQ] = {0};
    double dQ[nnzQ] = {0};

    for(int i = 0;i < curveNums;i++)
    {
        for (int j = 0; j < 21; j++)
        {
            if(j <= 5)
            {
                irowQ[j + 21*i] = j + 6 * i;   
            }
            else if(j > 5 && j < 11)
            {
                irowQ[j + 21*i] = (j - 5) + 6 * i;
            }
            else if(j > 10 && j < 15)
            {
                irowQ[j + 21*i] = (j - 9) + 6 * i;
            }
            else if(j > 14 && j < 18)
            {
                irowQ[j + 21*i] = (j - 12) + 6 * i;
            }
            else if(j > 17 && j < 20)
            {
                irowQ[j + 21*i] = (j - 14) + 6 * i;
            }
            else
            {
                irowQ[j + 21*i] = (j - 15) + 6 * i;
            }
        }
    }

    for(int i = 0;i < curveNums;i++)
    {
        for (int j = 0; j < 21; j++)
        {
            if(j <= 5)
            {
                jcolQ[j + 21*i] = irowQ[21*i];   
            }
            else if(j > 5 && j < 11)
            {
                jcolQ[j + 21*i] = irowQ[21*i + 1];
            }
            else if(j > 10 && j < 15)
            {
                jcolQ[j + 21*i] = irowQ[21*i + 2];
            }
            else if(j > 14 && j < 18)
            {
                jcolQ[j + 21*i] = irowQ[21*i + 3];
            }
            else if(j > 17 && j < 20)
            {
                jcolQ[j + 21*i] = irowQ[21*i + 4];
            }
            else
            {
                jcolQ[j + 21*i] = irowQ[21*i + 5];
            }
        }
    }

    for(int i = 0;i < curveNums;i++)
    {
        for (int j = 0; j < 21; j++)
        {
            if(j <= 5)
            {
                dQ[j + 21*i] = 0;   
            }
            else if(j > 5 && j < 11)
            {
                dQ[j + 21*i] = 0;
            }
            else if(j > 10 && j < 15)
            {
                dQ[j + 21*i] = 0;
            }
            else if(j > 14 && j < 18)
            {
                dQ[j + 21*i] = 0;
            }
            else if(j == 18)
            {
                dQ[j + 21*i] = 576 * std::pow(T(i),1);
            }
            else if(j == 19)
            {
                dQ[j + 21*i] = 1440 * std::pow(T(i),2);
            }
            else
            {
                dQ[j + 21*i] = 4800 * std::pow(T(i),3);
            }
        }
    }

    // for(int i = 0;i<nnzQ;i++)
    // {
    //     cout<<irowQ[i]<<" ";
    // }
    // cout<<endl;
    // for(int i=0;i<nnzQ;i++)
    // {
    //     cout<<jcolQ[i]<<" ";
    // }
    // cout<<endl;
    // for(int i = 0;i<nnzQ;i++)
    // {
    //     cout<<dQ[i]<<" ";
    // }
    // cout<<endl;

    int my = 4 * curveNums + 2;
    double b[my] = {0};
    for(int i = 0;i < my;i++)
    {
        if(i < 3)
        {
            b[i] = begin(i);
        }
        else if(i > 2 && i < (curveNums + 2))
        {
            b[i] = midlep(i - 3);
        }
        else if(i > (curveNums + 1) && i < (curveNums + 5))
        {
            b[i] = end(i - curveNums -2);
        }
        else
        {
            b[i] = 0;
        }
    }
    // for(int i = 0;i < my;i++)
    // {
    //     cout<<b[i]<<endl;
    // }
    
    //从这里开始测试所有数组的值是否正常。
    int nnzA = 36 * curveNums - 6;
    int irowA[nnzA] = {0};

    for(int i = 0;i < nnzA;)
    {
        int NonZeroRowNum = 3;
        int ContinueNum = 0;
        if(i < 6)
        {
            irowA[i] = 0;
            i++;
        }
        else if(i > 5 && i < 11)
        {
            irowA[i] = 1;
            i++;
        }
        else if(i > 10 && i < 15)
        {
            irowA[i] = 2;
            i++;
        }
        else if(i > 14 && i < (6 * curveNums + 8))
        {
            irowA[i] = NonZeroRowNum;
            irowA[i+1] = NonZeroRowNum;
            irowA[i+2] = NonZeroRowNum;
            irowA[i+3] = NonZeroRowNum;
            irowA[i+4] = NonZeroRowNum;
            irowA[i+5] = NonZeroRowNum;
            i = i + 6;
            NonZeroRowNum = NonZeroRowNum + 1;
        }
        else if(i > (6 * curveNums + 7) && i < (6 * curveNums + 14))
        {
            irowA[i] = NonZeroRowNum;
            irowA[i+1] = NonZeroRowNum;
            irowA[i+2] = NonZeroRowNum;
            irowA[i+3] = NonZeroRowNum;
            irowA[i+4] = NonZeroRowNum;
            irowA[i+5] = NonZeroRowNum;
            i = i+6;
            NonZeroRowNum = NonZeroRowNum + 1;
        }
        else if(i > (6 * curveNums + 13) && i < (6 * curveNums + 19))
        {
            irowA[i] = NonZeroRowNum;
            irowA[i+1] = NonZeroRowNum;
            irowA[i+2] = NonZeroRowNum;
            irowA[i+3] = NonZeroRowNum;
            irowA[i+4] = NonZeroRowNum;
            i = i + 5;
            NonZeroRowNum = NonZeroRowNum + 1;
        }
        else if(i > (6 * curveNums + 18) && i < (6 * curveNums + 23))
        {
            irowA[i] = NonZeroRowNum;
            irowA[i+1] = NonZeroRowNum;
            irowA[i+2] = NonZeroRowNum;
            irowA[i+3] = NonZeroRowNum;
            i = i + 4;
            NonZeroRowNum = NonZeroRowNum + 1;
        }

        else if(i > (6 * curveNums + 22 + ContinueNum * 30) && i < (6 * curveNums + 35 + ContinueNum * 30) && ContinueNum < (curveNums - 1))
        {
            for(int j = 0;j < 12;j++)
            {
                irowA[i + j] = NonZeroRowNum;
            }
            i = i + 12;
            NonZeroRowNum = NonZeroRowNum + 1;
        }

        else if(i > (6 * curveNums + 34 + ContinueNum * 30) && i < (6 * curveNums + 45 + ContinueNum * 30) && ContinueNum < (curveNums - 1))
        {
            for(int j = 0;j < 10;j++)
            {
                irowA[i + j] = NonZeroRowNum;
            }
            i = i + 10;
            NonZeroRowNum = NonZeroRowNum + 1;
        }

        else if(i > (6 * curveNums + 44 + ContinueNum * 30) && i < (6 * curveNums + 53 + ContinueNum * 30) && ContinueNum < (curveNums - 1))
        {
            for(int j = 0;j < 8;j++)
            {
                irowA[i + j] = NonZeroRowNum;
            }
            i = i + 8;
            NonZeroRowNum = NonZeroRowNum + 1;
            ContinueNum = ContinueNum + 1;
        }
    }

    int jcolA[nnzA] = {0};
    for(int i = 0;i < nnzA;)
    {
        int PNum = 0;
        int NonContinueNum = 0;
        int ContinueNum = 0;
        int NonZeroColNum = 0;
        if(i < 6)
        {
            for(int j = 0;j < 6;j++)
            {
                jcolA[i + j] = NonZeroColNum + j;
            }
            i = i + 6;
        }
        else if(i > 5 && i < 11)
        {
            for(int j = 0;j < 5;j++)
            {
                jcolA[i + j] = NonZeroColNum + j + 1;
            }
            i = i + 5;
        }
        else if(i > 10 && i < 15)
        {
            for(int j = 0;j < 4;j++)
            {
                jcolA[i + j] = NonZeroColNum + j + 2;
            }
            i = i + 4;
        }
        else if(i > (14 + 6 * PNum) && i < (21 + 6 * PNum) && PNum < (curveNums - 1))
        {
            for(int j = 0;j < 6;j++)
            {
                jcolA[i + j] = NonZeroColNum + j + 6 * (PNum + 1);
            }
            i = i + 6;
            PNum = PNum + 1;
        }
        //15 + 6 * (k - 1)
        else if(i > (6 * curveNums + 8) && i < (6 * curveNums + 15))
        {
            for(int j = 0;j < 6;j++)
            {
                jcolA[i + j] = NonZeroColNum + j + 6 * curveNums;
            }
            i = i + 6;
        }
        //i = 37
        else if(i > (6 * curveNums + 14) && i < (6 * curveNums + 20))
        {
            for(int j = 0;j < 5;j++)
            {
                jcolA[i + j] = NonZeroColNum + j + 1 + 6 * curveNums;
            }
            i = i + 5;
        }
        //i = 42
        else if(i > (6 * curveNums + 19) && i < (6 * curveNums + 24))
        {
            for(int j = 0;j < 4;j++)
            {
                jcolA[i + j] = NonZeroColNum + j + 2 + 6 * curveNums;
            }
            i = i + 4;
        }

        else if(i > (6 * curveNums + 23 + ContinueNum * 30) && i < (6 * curveNums + 36 + ContinueNum * 30) && ContinueNum < (curveNums - 1))
        {
            for(int j = 0;j < 12;j++)
            {
                jcolA[i + j] = NonZeroColNum + j + 6 * (ContinueNum + 1);
            }
            i = i + 12;
        }

        else if(i > (6 * curveNums + 35 + ContinueNum * 30) && i < (6 * curveNums + 46 + ContinueNum * 30) && ContinueNum < (curveNums - 1))
        {
            for(int j = 0;j < 5;j++)
            {
                jcolA[i + j] = NonZeroColNum + j + 1 + 6 * (ContinueNum + 1); 
            }
            for(int j = 5;j < 10;j++)
            {
                jcolA[i + j] = NonZeroColNum + j + 1 + 6 * (ContinueNum + 1);
            }
            i = i + 10;
        }

        else if(i > (6 * curveNums + 45 + ContinueNum * 30) && i < (6 * curveNums + 54 + ContinueNum * 30) && ContinueNum < (curveNums - 1))
        {
            for(int j = 0;j < 4;j++)
            {
                jcolA[i + j] = NonZeroColNum + j + 2 + 6 * (ContinueNum + 1); 
            }
            for(int j = 4;j < 8;j++)
            {
                jcolA[i + j] = NonZeroColNum + j + 2 + 6 * (ContinueNum + 1);
            }
            i = i + 8;
            ContinueNum = ContinueNum + 1;
        }
    }

    double dA[nnzA] = {0};
    for(int i = 0;i < nnzA;)
    {
        int ThCurve = 0;
        int ThConCurve = 0;
        int PNum = 0;
        int NonContinueNum = 0;
        int ContinueNum = 0;
        int NonZeroColNum = 0;
        if(i < 6)
        {
            for(int j = 0;j < 6;j++)
            {
                jcolA[i + j] = std::pow(T(ThCurve),j);
            }
            i = i + 6;
        }
        else if(i > 5 && i < 11)
        {
            for(int j = 0;j < 5;j++)
            {
                jcolA[i + j] = (j + 1) * std::pow(T(ThCurve),j);
            }
            i = i + 5;
        }
        else if(i > 10 && i < 15)
        {
            for(int j = 0;j < 4;j++)
            {
                jcolA[i + j] = (j + 2) * std::pow(T(ThCurve),j);
            }
            i = i + 4;
            ThCurve = ThCurve + 1;
        }
        else if(i > (14 + 6 * PNum) && i < (21 + 6 * PNum) && PNum < (curveNums - 1))
        {
            for(int j = 0;j < 6;j++)
            {
                jcolA[i + j] = std::pow(T(ThCurve),j);
            }
            i = i + 6;
            PNum = PNum + 1;
            ThCurve = ThCurve + 1;
        }
        //15 + 6 * (k - 1)
        else if(i > (6 * curveNums + 8) && i < (6 * curveNums + 15))
        {
            for(int j = 0;j < 6;j++)
            {
                jcolA[i + j] = std::pow(T(ThCurve),j);
            }
            i = i + 6;
        }
        //i = 37
        else if(i > (6 * curveNums + 14) && i < (6 * curveNums + 20))
        {
            for(int j = 0;j < 5;j++)
            {
                jcolA[i + j] = (j + 1) * std::pow(T(ThCurve),j);
            }
            i = i + 5;
        }
        //i = 42
        else if(i > (6 * curveNums + 19) && i < (6 * curveNums + 24))
        {
            for(int j = 0;j < 4;j++)
            {
                jcolA[i + j] = (j + 2) * std::pow(T(ThCurve),j);
            }
            i = i + 4;
            ThCurve = ThCurve + 1;
        }

        else if(i > (6 * curveNums + 23 + ContinueNum * 30) && i < (6 * curveNums + 36 + ContinueNum * 30) && ContinueNum < (curveNums - 1))
        {
            for(int j = 0;j < 12;j++)
            {
                if(j < 6)
                {
                    jcolA[i + j] = std::pow(T(ThConCurve),j);
                }
                else
                {
                    if(j == 6)
                    {
                        jcolA[i + j] = -1;
                    }
                    else
                    {
                        // jcolA[i + j] = -0 * std::pow(T(ThConCurve),(j - 6));
                        jcolA[i + j] = 0;
                    }
                }
            }
            i = i + 12;
        }

        else if(i > (6 * curveNums + 35 + ContinueNum * 30) && i < (6 * curveNums + 46 + ContinueNum * 30) && ContinueNum < (curveNums - 1))
        {
            for(int j = 0;j < 5;j++)
            {
                jcolA[i + j] = (j + 1) * std::pow(T(ThConCurve),j);
            }
            for(int j = 5;j < 10;j++)
            {
                if(j == 5)
                {
                    jcolA[i + j] = -1;
                }
                else
                {
                    // jcolA[i + j] = -0;
                    jcolA[i + j] = 0;
                }
            }
            i = i + 10;
        }

        else if(i > (6 * curveNums + 45 + ContinueNum * 30) && i < (6 * curveNums + 54 + ContinueNum * 30) && ContinueNum < (curveNums - 1))
        {
            for(int j = 0;j < 4;j++)
            {
                jcolA[i + j] = (j + 2) * std::pow(T(ThConCurve),j); 
            }
            for(int j = 4;j < 8;j++)
            {
                if(j == 4)
                {
                    jcolA[i + j] = -2;
                }
                else
                {
                    jcolA[i + j] = 0;
                }
            }
            i = i + 8;
            ContinueNum = ContinueNum + 1;
            ThConCurve = ThConCurve + 1;
        }
    }

    const int mz = 0;
    double *clow = 0;
    char *iclow = 0;
    double *cupp = 0;
    char *icupp = 0;

    const int nnzC = 0;
    int *irowC = 0;
    int *jcolC = 0;
    double *dC = 0;

    QpGenSparseMa27 * qp = new QpGenSparseMa27(nx,my,mz,nnzQ,nnzA,nnzC);
    QpGenData *prob = (QpGenData *)qp ->copyDataFromSparseTriple(
                                                                                                                                                c,irowQ,nnzQ,jcolQ,dQ,
                                                                                                                                                xlow,ixlow,xupp,ixupp,
                                                                                                                                                irowA,nnzA,jcolA,dA,b,
                                                                                                                                                irowC,nnzC,jcolC,dC,
                                                                                                                                                clow,iclow,cupp,icupp);
    QpGenVars * vars = (QpGenVars *)qp ->makeVariables(prob);
    QpGenResiduals *resid = (QpGenResiduals *)qp -> makeResiduals(prob);
    GondzioSolver *s = new GondzioSolver(qp,prob);

    s->monitorSelf();
    int ierr = s -> solve(prob,vars,resid);
    if(ierr == 0)
    {
        cout.precision(4);
        cout<<"Solution: "<<endl;
        vars->x->writefToStream(cout,"x[%{index}] = %{value}");
    }
    else
    {
        cout<<"Could not solve the problem !"<<endl;
    }
}
