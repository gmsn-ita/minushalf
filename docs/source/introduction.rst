=============
Introduction
=============


Uma visão intuitiva do método DFT -1/2
-------------------------------------------
DFT-1/2, an alternative way of referring to the LDA -1/2 [1]_ [2]_ and GGA -1/2 [2]_ techniques, 
is a method that performs semiconductor band-gap calculations with precision close 
to the state of the art algorithms [2]_. These technique aim to expand the half-occupation 
technique [3]_ [4]_ [5]_, formalized by Janak's theorem, to crystals using modern exchange-correlation approaches [6]_ [7]_.

A intuição do método vem do fato que as bandas de energia de um cristal são formadas
pelo overlap de orbitais atômicos, principalmente aqueles que  formam as camadas mais externas [8]_. Essa
influência pode ser quantificada pela projeção da função de  onda em um dado orbital, a Figura 1 mostra
o caráter da última banda de valência e da primeira banda de condução para o CdO, a cor magneta representa
o o caráter para o tipo de orbital d e o amarelo o caráter do orbital p. Assim, uma meia ocupação dos autovalores
a nível atômico poderia se propagar e resultar em uma correção no Gap para sistemas cristalinos.

.. figure:: images/cdo_bands.png
   :width: 500

   Caráter do orbital para as bandas de valência do CdO. O caráter p é representado
   em amarelo e o caráter d em magneta.


Diferentemente dos átomos, onde os níveis de energia são discretos, buscou-se aplicar a meia ocupação
nas bandas do cristal. Para isso, não pode-se fazer a meia ocupação diretamente pela mudança de 
densidade, uma vez que seria gerado um sistema infinitamente carregado, o que levaria a uma 
divergência nos cálculos. Também seria irrelevante conseguir modificar somente uma quantidade finita
de elétrons em células unitárias, uma vez que a carga iria se tornar irrelevante para a quantidade infinita de 
células unitárias. Para contornar esse problema, são feitas pequenas aproximações onde se torna possível
encontrar o potencial meio ocupado do cristal através de outros potenciais, assim como mostra a equação
abaixo;

.. math:: 
   V_{crystal}^{-1/2} = V_{crystal} - V_{1/2e}

Onde :math:`V_{crystal}^{-1/2}` é o potencial do cristal meio ocupado, :math:`V_{crystal}`
é o potencial do cristal com a ocupação padrão e :math:`V_{1/2e}` é o potencial do respectivo nível
ocupado com meio elétron.

Para gerar :math:`V_{1/2e}`, utilizamos a seguinte equação para os respectivos átomos do cristal

.. math:: 
   V_{1/2e} = V_{atom} - V_{atom}^{f_{\alpha}=-1/2}


Onde :math:`V_{atom}` é o potencial do átomo com a ocupação padrão, :math:`V_{atom}^{f_{\alpha}=-1/2}`
é o potencial do átomo com o nível :math:`\alpha` ocupado.

Como os átomos se repetem o em cada célula primitiva, o potencial :math:`V_{1/2e}` é periódico, somando
com o fato de que :math:`V_{crystal}` é periódico, pode-se concluir que :math:`V_{crystal}^{-1/2}` 
é periódico.

Porém, ainda resta um problema. Como é gerado um sistema artificial de carga não nula, a correção
em uma célula gera um potencial artificial que deveria ser somado nas células vizinhas, o que levaria
a uma divergência no potencial, uma vez que temos infinitas células unitárias. Para essa finalidade,
adicionou-se um parâmetro de corte para esse potencial conforme a equação abaixo

.. math::
   V_{1/2e} = (V_{atom} - V_{atom}^{f_{\alpha}=-1/2})\cdot \theta (r)

.. math::
   \theta (r) = \left\{\begin{matrix}
   \theta (r) = A \cdot[1-(\frac{r}{CUT})^{8}]^{3}) , r \leq CUT \\
   \theta (r) = 0 , r > CUT
   \end{matrix}\right.

Onde CUT é o raio de corte e A é um fator de escala chamado amplitude. Por meio de argumentos variacionais
,o cut e a amplitude são escolhidos a partir da maximização do Gap.

Onde realizar a meia ocupação?
--------------------------------------
A meia ocupação deve ser realizada na última banda de valência e na primeira banda da condução. Existem dois tipos de 
correção, a simples e a fracionária. A escolha de qual correção fazer depende da análise da composição da banda. Suponha que
temos uma matriz onde os átomos da célula unitária são representados como linhas e os tipos de orbitais atômicos como coluna,
cada valor `a_{ij}` representa, em porcentagem, quanto aquele dado orbital de um determinado átomo contribui para a formação da
banda.

.. math ::
   A = \begin{bmatrix} 
   a_{11} & a_{12} & \dots \\
   \vdots & \ddots & \\
   a_{N1} &        & a_{NK} 
   \end{bmatrix}

Onde:

.. math ::
   \sum_{i=1}^{N} \sum_{j=1}^{K} a_{ij} = 100



Correção simples
#########################
Aplica-se o método da correção simples quando um índice :math:`a_{ij}` representa majoritariamente a 
composição da banda, de forma que a influência dos demais orbitais é desprezível. 
Assim, a correçao de meio elétron é feita somente no orbital :math:`j` do átomo :math:`i`.

Correção fracionária
########################
Aplica-se o método da correção fracionária quando diferentes orbitais atômicos tem influência siguinificativa
na composição da banda. Para distribuir o meio elétron, escolhe-se um treshold
:math:`\epsilon` , que representa o valor mínimo de :math:`a_{ij}` considerado na correção. Definidos esses
valores, o meio elétron será dividido entre os átomos, de forma proporcional ao coeficiente :math:`a_{ij}`.

Considerações finais
#######################

Após aplicada a correção, deve-se encontrar o cut e amplitude ótimos para cada átomo corrigido para, finalmente,
encontrarmos o valor final para o Gap.

Vale salientar que , em muitos casos, a correção na banda de valência já retorna resultados suficientes e próximos o suficiente
dos valores experimental, o que descarta a necessidade de uma correção adicional na banda de condução.


DFT -1/2 results
----------------------

The results obtained by the application the method has the same precision [2]_ as the GW [9]_ algorithm , considered 
the state of the art for calculating the band-gap of semiconductors. In addition, the computational complexity of the method 
is equivalent to calculating the Khon-Shan gap, which allows the technique to be applied to complex systems.


References
------------------

.. [1] L. G. Ferreira, M. Marques, and L. K. Teles, `Phys. Rev. B 78, 125116 (2008) <http://dx.doi.org/10.1103/PhysRevB.78.125116>`_.

.. [2] L. G. Ferreira, M. Marques, and L. K. Teles, `AIP Adv. 1, 032119 (2011) <https://doi.org/10.1063/1.3624562>`_.

.. [3] J.C. Slater and K. H. Johnson, `Phys. Rev. B 5, 844 (1972) <http://dx.doi.org/10.1103/PhysRevB.5.844>`_.

.. [4] J.C. Slater, `Adv. Quantum Chem. 6, 1 (1972) <http://dx.doi.org/10.1016/S0065-3276(08)60541-9>`_.

.. [5] J. C. Slater and J. H. Wood, Int. J. Quant. Chem. Suppl. 4, 3 (1971).

.. [6] J. P. Perdew and A. Zunger, `Phys. Rev. B 23, 5048 (1981) <http://dx.doi.org/10.1103/PhysRevB.23.5048>`_.

.. [7] J. P. Perdew, K. Burke, and M. Ernzerhof, `Phys. Rev. Lett. 77, 3865 (1996) <http://dx.doi.org/10.1103/PhysRevLett.77.3865>`_ .

.. [8] Holgate, Sharon Ann (2009). Understanding Solid State Physics. CRC Press. pp. 177–178. ISBN 978-1-4200-1232-3.

.. [9] G. Onida, L. Reining, and A. Rubio, `Rev. Mod. Phys. 74, 601 (2002) <http://dx.doi.org/10.1103/RevModPhys.74.601>`_.
