import numpy as np
import matplotlib.pyplot as plt 
import sys
epss=10e-15
print("晶体学可视化程序 Release v1.0")
print("Edited by Python")
print("Produced by Chemyy")
while(1):
    print("等待输入")
    print("1:主程序")
    print("2:说明文件")
    print("3:退出")
    inp=input("输入命令(三选一):")
    if inp=="2":
        print("使用说明书")
        print("本程序基于python的numpy库实现矩阵规范化运算、基于matplotlib库实现晶胞可视化,并移植到windows的cmd平台,可直接以exe的格式打开。")
        print("目前能够实现的功能有六项,分别是1.绘制平移格子 2.绘制原胞和基向量 3.绘制W-S晶胞的顶点 4.绘制倒易平移格子 5.绘制倒易原胞和倒易基矢量 6.绘制W-S晶胞和第一布里渊区")
        print("其中4、6两项利用解线性方程组的方法获得顶点,由于numpy运算精度有限,可能存在误判现象,这时候建议将基矢量长度增加。")
        print("关于基矢量的输入,此处采用原胞输入而非惯用正当晶胞,现给出14种布拉维格子的原胞关系:")
        print("1.三斜:a!=b!=c,α!=β!=γ")
        print("2.单斜P:a!=b!=c,α=γ=90")
        print("3.单斜C:a=b!=c,α!=β!=γ")
        print("4.正交P:a!=b!=c,α=β=γ=90")
        print("5.正交C :a=b!=c,α=β=90")
        print("6.正交I:a=b=c,α!=β!=γ")
        print("7.正交F:a!=b!=c,a=bcos(γ)+ccos(β)")
        print("8.六方P:a=b,α=β=90,γ=120")
        print("9.三方R:a=b=c,α=β=γ")
        print("10.四方P:a=b!=c,α=β=γ=90")
        print("11.四方I:a=b=c,α=β!=γ")
        print("12.立方P:a=b=c,α=β=γ=90")
        print("13.立方I,a=b=c,α=β=γ=109.5")
        print("14.立方F:a=b=c,α=β=γ=60")
        print("在可视化界面,按住左键移动实现拖动视角, 按住右键上下移动实现缩放,按住中键移动实现平移")
        print("")
        print("更新日志")
        print("v202309272349:实现格点平移、基矢量、原胞的绘制")
        print("v202310071636:实现W-S晶胞和第一布里渊区顶点的绘制")
        print("v202310072042:修了矩阵奇异的bug,修缮了说明文件")
        print("v202310092017:给WSC和FBZ加了线框,提升了用户交互体验")
        print("")
        print("已知问题")
        print("W-S晶胞和第一布里渊区无法预测地可能发生缺线或多线的情况")
        input("输入任意键返回...")
        continue
    if inp=="3":
        sys.exit()
    if inp=="1":
        ax = plt.axes(111,projection='3d', proj_type='ortho')
        a=float(input("原胞第一基矢量比例长度(无量纲,建议大于等于1小于等于10)a="))
        b=float(input("原胞第二基矢量比例长度(无量纲,建议大于等于1小于等于10)b="))
        c=float(input("原胞第三基矢量比例长度(无量纲,建议大于等于1小于等于10)c="))
        alpha=float(input("第一第二基矢量交角(单位:度)alpha="))/180*np.pi
        beta=float(input("第二第三基矢量交角(单位:度)bata="))/180*np.pi
        gamma=float(input("第三第一基矢量交角(单位:度)gamma="))/180*np.pi
        a_1=np.array([a,0,0])                                                     
        a_2=np.array([b*np.cos(gamma),b*np.sin(gamma),0])
        a_3=np.array([c*np.cos(beta),c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma),c*np.sqrt(1-np.cos(beta)**2-((np.cos(alpha)-np.cos(beta)*np.cos(gamma))/np.sin(gamma))**2)])
        Real_Bases=np.vstack((np.vstack((a_1,a_2)),a_3))
        V_cell=np.linalg.det(Real_Bases)
        if V_cell==0:
            print("错误:实空间基底共面,请检查输入。")
            continue
        else:
            b_1=2*np.pi/V_cell*np.cross(a_2,a_3)
            b_2=2*np.pi/V_cell*np.cross(a_3,a_1)
            b_3=2*np.pi/V_cell*np.cross(a_1,a_2)
        Reciprocal_Bases=np.vstack((np.vstack((b_1,b_2)),b_3))
        ita=1.5
        ax.set_xlim3d([-ita*a,ita*a])
        ax.set_ylim3d([-ita*b,ita*b])
        ax.set_zlim3d([-ita*c,ita*c])
        ax.set_aspect('equal')
        ax.set_axis_off()
        print("等待绘图指令……")
        print("1:实空间点阵(蓝色格点)")
        print("2:实空间基向量及原胞格子(蓝色向量与线框)")
        print("3:W-S晶胞(绿色线框)")
        print("4:倒空间点阵(红色格点)")
        print("5:倒空间基向量及原胞格子(红色向量与线框)")
        print("6:第一布里渊区(黄色线框)")
        order=input("空格隔开需求对应的序号(建议一次至多三个):").split()
        if "1" in order:
            Scale=int((int(input("实点阵左右拓展个数(奇数,建议取1):"))+1)/2)
            for i in range(-Scale,Scale+1):
                for j in range(-Scale,Scale+1):
                    for k in range(-Scale,Scale+1):
                        vr=i*a_1+j*a_2+k*a_3
                        ax.plot3D(vr[0],vr[1],vr[2],c='blue',marker='o')
        if "2" in order:
            ax.quiver(0,0,0,Real_Bases[0][0],Real_Bases[0][1],Real_Bases[0][2],color="blue") 
            ax.quiver(0,0,0,Real_Bases[1][0],Real_Bases[1][1],Real_Bases[1][2],color="blue")
            ax.quiver(0,0,0,Real_Bases[2][0],Real_Bases[2][1],Real_Bases[2][2],color="blue")
            ax.text(Real_Bases[0][0],Real_Bases[0][1],Real_Bases[0][2]+0.1,"$a_1$",color="blue") 
            ax.text(Real_Bases[1][0],Real_Bases[1][1],Real_Bases[1][2]+0.1,"$a_2$",color="blue")
            ax.text(Real_Bases[2][0],Real_Bases[2][1],Real_Bases[2][2]+0.1,"$a_3$",color="blue")
            dx=[]
            dy=[]
            dz=[]
            for i in [0,1]:
                for j in [0,1]:
                    for k in[0,1]:
                        vs=i*a_1+j*a_2+k*a_3
                        dx.append(vs[0])
                        dy.append(vs[1])
                        dz.append(vs[2])
            A,B,C,D,E,F,G,H=zip(dx,dy,dz)
            L= zip(A,E,G,C,A,B,F,H,D,B,F,E,F,H,G,H,D,C,D)
            ax.plot3D(*L,zdir='z',c='b',ms=10,alpha=0.5)
        if "3" in order:
            Normal_Vector=np.array([0,0,0])
            Vertex=np.array([[0,0,0],[0,0,0]])
            for i in [-1,0,1]:
                for j in [-1,0,1]:
                    for k in [-1,0,1]:
                        vr=i*a_1+j*a_2+k*a_3
                        Normal_Vector=np.vstack((Normal_Vector,vr))
            for i in range(1,Normal_Vector.shape[0]):
                for j in range(i+1,Normal_Vector.shape[0]):
                    for k in range(j+1,Normal_Vector.shape[0]):
                        M=Normal_Vector[i]
                        M=np.vstack((M,Normal_Vector[j]))
                        M=np.vstack((M,Normal_Vector[k]))
                        B=(np.linalg.norm(M[0]))**2
                        B=np.vstack((B,(np.linalg.norm(M[1]))**2))
                        B=np.vstack((B,(np.linalg.norm(M[2]))**2))
                        B=0.5*B
                        if np.linalg.det(M)!=0:
                            M=np.linalg.inv(M)
                            r=M.dot(B)
                        else:
                            continue
                        delta=1
                        for p in range(1,Normal_Vector.shape[0]):
                            RU=r[0]*Normal_Vector[p][0]+r[1]*Normal_Vector[p][1]+r[2]*Normal_Vector[p][2]
                            U2=(Normal_Vector[p][0])**2+Normal_Vector[p][1]**2+Normal_Vector[p][2]**2
                            if RU[0]-0.5*U2>epss:
                                delta=0
                        if delta==1:
                            dell=0
                            for p in range(1,Vertex.shape[0]):
                                if (abs(r[0]-Vertex[p][0])<epss)&(abs(r[1]-Vertex[p][1])<epss)&(abs(r[2]-Vertex[p][2])<epss):
                                    dell+=1
                            if dell==0:
                                Vertex=np.vstack((Vertex,r.T))
            L=np.array([0])
            for i in range(1,Normal_Vector.shape[0]):
                for j in range(i+1,Normal_Vector.shape[0]):
                    M=Normal_Vector[i]
                    M=np.vstack((M,Normal_Vector[j]))
                    tmp0=Normal_Vector[i][1]*Normal_Vector[j][2]-Normal_Vector[i][2]*Normal_Vector[j][1]
                    tmp1=Normal_Vector[i][2]*Normal_Vector[j][0]-Normal_Vector[i][0]*Normal_Vector[j][2]
                    tmp2=Normal_Vector[i][0]*Normal_Vector[j][1]-Normal_Vector[i][1]*Normal_Vector[j][0]
                    tmp=0.25*np.array([tmp0,tmp1,tmp2])
                    M=np.vstack((M,tmp))
                    sigma0=0.5*np.linalg.norm(Normal_Vector[i])**2
                    sigma1=0.5*np.linalg.norm(Normal_Vector[j])**2
                    B=np.array([[sigma0],[sigma1],[0]])
                    if np.linalg.det(M)!=0:
                        M=np.linalg.inv(M)
                        r=M.dot(B)
                        l=np.linalg.norm(r)
                        L=np.vstack((L,l))
            for i in range(2,Vertex.shape[0]):
                ax.scatter3D(Vertex[i][0],Vertex[i][1],Vertex[i][2],c='green',marker='v')
            for i in range(2,Vertex.shape[0]):
                for j in range(i+1,Vertex.shape[0]):
                    x1=Vertex[i][0]
                    x2=Vertex[j][0]
                    y1=Vertex[i][1]
                    y2=Vertex[j][1]
                    z1=Vertex[i][2]
                    z2=Vertex[j][2]
                    ldb=1.0*(x1*(x1-x2)+y1*(y1-y2)+z1*(z1-z2))/((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
                    d=np.array([x1,y1,z1])+ldb*np.array([x2-x1,y2-y1,z2-z1])
                    dl=np.linalg.norm(d)
                    for p in range(1,L.shape[0]):
                        if abs(dl-L[p])<=epss:
                            ax.plot3D([x1,x2],[y1,y2],[z1,z2],zdir='z',c='green',ms=10,alpha=0.5)
                            break
        if "4" in order:
            Scale=int((int(input("倒点阵左右拓展个数(奇数,建议取1):"))+1)/2)
            for i in range(-Scale,Scale+1):
                for j in range(-Scale,Scale+1):
                    for k in range(-Scale,Scale+1):
                        vk=i*b_1+j*b_2+k*b_3
                        ax.plot3D(vk[0],vk[1],vk[2],c='red',marker='*')
        if "5" in order:
            ax.quiver(0,0,0,Reciprocal_Bases[0][0],Reciprocal_Bases[0][1],Reciprocal_Bases[0][2],color="red") 
            ax.quiver(0,0,0,Reciprocal_Bases[1][0],Reciprocal_Bases[1][1],Reciprocal_Bases[1][2],color="red")
            ax.quiver(0,0,0,Reciprocal_Bases[2][0],Reciprocal_Bases[2][1],Reciprocal_Bases[2][2],color="red")
            ax.text(Reciprocal_Bases[0][0],Reciprocal_Bases[0][1],Reciprocal_Bases[0][2]+0.1,"$b_1$",color="red") 
            ax.text(Reciprocal_Bases[1][0],Reciprocal_Bases[1][1],Reciprocal_Bases[1][2]+0.1,"$b_2$",color="red")
            ax.text(Reciprocal_Bases[2][0],Reciprocal_Bases[2][1],Reciprocal_Bases[2][2]+0.1,"$b_3$",color="red")
            dx=[]
            dy=[]
            dz=[]
            for i in [0,1]:
                for j in [0,1]:
                    for k in[0,1]:
                        vs=i*b_1+j*b_2+k*b_3
                        dx.append(vs[0])
                        dy.append(vs[1])
                        dz.append(vs[2])
            A,B,C,D,E,F,G,H=zip(dx,dy,dz)
            L= zip(A,E,G,C,A,B,F,H,D,B,F,E,F,H,G,H,D,C,D)
            ax.plot3D(*L,zdir='z',c='r',ms=10,alpha=0.5)
        if "6" in order:
            Normal_Vector=np.array([0,0,0])
            Vertex=np.array([[0,0,0],[0,0,0]])
            for i in [-1,0,1]:
                for j in [-1,0,1]:
                    for k in [-1,0,1]:
                        vr=i*b_1+j*b_2+k*b_3
                        Normal_Vector=np.vstack((Normal_Vector,vr))
            for i in range(1,Normal_Vector.shape[0]):
                for j in range(i+1,Normal_Vector.shape[0]):
                    for k in range(j+1,Normal_Vector.shape[0]):
                        M=Normal_Vector[i]
                        M=np.vstack((M,Normal_Vector[j]))
                        M=np.vstack((M,Normal_Vector[k]))
                        B=(np.linalg.norm(M[0]))**2
                        B=np.vstack((B,(np.linalg.norm(M[1]))**2))
                        B=np.vstack((B,(np.linalg.norm(M[2]))**2))
                        B=0.5*B
                        if np.linalg.det(M)!=0:
                            M=np.linalg.inv(M)
                            r=M.dot(B)
                        else:
                            continue
                        delta=1
                        for p in range(1,Normal_Vector.shape[0]):
                            RU=r[0]*Normal_Vector[p][0]+r[1]*Normal_Vector[p][1]+r[2]*Normal_Vector[p][2]
                            U2=(Normal_Vector[p][0])**2+Normal_Vector[p][1]**2+Normal_Vector[p][2]**2
                            if RU[0]-0.5*U2>epss:
                                delta=0
                        if delta==1:
                            dell=0
                            for p in range(1,Vertex.shape[0]):
                                if (abs(r[0]-Vertex[p][0])<epss)&(abs(r[1]-Vertex[p][1])<epss)&(abs(r[2]-Vertex[p][2])<epss):
                                    dell+=1
                            if dell==0:
                                Vertex=np.vstack((Vertex,r.T))
            L=np.array([0])
            for i in range(1,Normal_Vector.shape[0]):
                for j in range(i+1,Normal_Vector.shape[0]):
                    M=Normal_Vector[i]
                    M=np.vstack((M,Normal_Vector[j]))
                    tmp0=Normal_Vector[i][1]*Normal_Vector[j][2]-Normal_Vector[i][2]*Normal_Vector[j][1]
                    tmp1=Normal_Vector[i][2]*Normal_Vector[j][0]-Normal_Vector[i][0]*Normal_Vector[j][2]
                    tmp2=Normal_Vector[i][0]*Normal_Vector[j][1]-Normal_Vector[i][1]*Normal_Vector[j][0]
                    tmp=0.25*np.array([tmp0,tmp1,tmp2])
                    M=np.vstack((M,tmp))
                    sigma0=0.5*np.linalg.norm(Normal_Vector[i])**2
                    sigma1=0.5*np.linalg.norm(Normal_Vector[j])**2
                    B=np.array([[sigma0],[sigma1],[0]])
                    if np.linalg.det(M)!=0:
                        M=np.linalg.inv(M)
                        r=M.dot(B)
                        l=np.linalg.norm(r)
                        L=np.vstack((L,l))
            for i in range(2,Vertex.shape[0]):
                ax.scatter3D(Vertex[i][0],Vertex[i][1],Vertex[i][2],c='brown',marker='x')
            for i in range(2,Vertex.shape[0]):
                for j in range(i+1,Vertex.shape[0]):
                    x1=Vertex[i][0]
                    x2=Vertex[j][0]
                    y1=Vertex[i][1]
                    y2=Vertex[j][1]
                    z1=Vertex[i][2]
                    z2=Vertex[j][2]
                    ldb=1.0*(x1*(x1-x2)+y1*(y1-y2)+z1*(z1-z2))/((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
                    d=np.array([x1,y1,z1])+ldb*np.array([x2-x1,y2-y1,z2-z1])
                    dl=np.linalg.norm(d)
                    for p in range(1,L.shape[0]):
                        if abs(dl-L[p])<=epss:
                            ax.plot3D([x1,x2],[y1,y2],[z1,z2],zdir='z',c='brown',ms=10,alpha=0.5)
                            break
        print("少女祈祷中...")
        plt.show()
    else:
        print("请检查你的输入")
        continue