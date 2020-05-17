#include <vector>
#include <fstream>
#include <stdio.h>
#include <cstdio>
#include <set>
#include <map>
#include <algorithm>
#include <QTextStream>
#include <QDebug>
#include <QFile>
#include <eigen/Eigen/Core>
#include <eigen/Eigen/IterativeLinearSolvers>
#include "placer.h"


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> DMat;

struct SubRegion
{
    double top = 100.0;
    double bot = 0.0;
    double left = 0.0;
    double right = 100.0;
    std::shared_ptr<SubRegion> sr1;
    std::shared_ptr<SubRegion> sr2;
    GateVector gates;
    GateVector pads;
    bool placed = false;

    double height() const { assert(top > bot); return top - bot; }
    double width() const { assert(right > left); return right - left; }
};

struct Edge 
{
    int id1 = -1; // id of first gate or pad to connect to
    int id2 = -1; // id of first gate or pad to connect to
    double weight = 0.0; // weight of this connection
};

struct Gate
{
    int id = -1;
    double x = 0.0;
    double y = 0.0;
    std::vector<int> nets;
};

std::vector<int> splitLine(const std::string& line)
{
    std::vector<int> output;
    std::string word = "";
    for(char c : line)
    {
        if(isspace(c))
        {
            if(word.size() > 0)
            {
                output.push_back(atoi(word.c_str()));
            }
            word = "";
        }
        else
        {
            word.push_back(c);
        }
    }
    if(word.size() > 0)
    {
        output.push_back(atoi(word.c_str()));
    }
    return output;
}

std::shared_ptr<Gate> get(const GateVector& gates, int id)
{
    for(size_t ii = 0; ii < gates.size(); ii++)
    {
        if(gates[ii]->id == id)
        {
            return gates[ii];
        }
    }
    return nullptr;
}

bool readNetlist(const char* filename, GateVector& gates, GateVector& pads)
{
    QFile input(filename);
    if(!input.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        printf("cannot open file '%s'\n", filename);
        return false;
    }

    std::vector<std::vector<int>> lines;
    while(!input.atEnd())
    {
        lines.push_back(splitLine(input.readLine().toStdString()));
    }
    if(lines.size() < 1)
    {
        printf("empty input\n");
        return false;
    }
    if(lines[0].size() != 2)
    {
        printf("first line: expected: [numGates] [numNets]\n");
        return false;
    }
    int pos = 0;
    int numGates = lines[pos][0];
    int numNets = lines[pos][1];
    pos++;

    if((lines.size() - pos) < size_t(numGates))
    {
        printf("malformatted input: not enough gates found\n");
        return false;
    }

    for(int ii = 0; ii < numGates; ii++)
    {
        const std::vector<int>& thisLine = lines[pos++];

        if(thisLine.size() < 2)
        {
            printf("line %d malformatted\n", pos);
            return false;
        }

        int connections = thisLine[1];
        if(thisLine.size() != size_t(2 + connections))
        {
            printf("line %d malformatted\n", pos);
            return false;
        }

        auto g = std::make_shared<Gate>();
        g->id = thisLine[0];
        g->nets.resize(connections);
        for(int jj = 0; jj < connections; jj++)
        {
           g->nets[jj] = thisLine[2 + jj];
        }
        gates.push_back(g);
    }

    if((size_t(pos) >= lines.size()) || (lines[pos].size() != 1))
    {
        printf("malformatted input: expected num pads at line %d\n", pos);
        return false;
    }
    int numPads = lines[pos++][0];

    if((lines.size() - pos) < size_t(numPads))
    {
        printf("malformatted input: not enough pads found\n");
        return false;
    }

    for(int ii = 0; ii < numPads; ii++)
    {
        const std::vector<int>& thisLine = lines[pos++];

        if(thisLine.size() != 4)
        {
            printf("line %d malformatted\n", pos);
            return false;
        }

        auto g = std::make_shared<Gate>();
        g->id = thisLine[0] + numGates;
        g->nets.push_back(thisLine[1]);
        g->x = thisLine[2];
        g->y = thisLine[3];
        pads.push_back(g);
    }
    return true;
}

void writeNetlist(const char* filename, const GateVector& gates, const GateVector& pads)
{
    std::set<int> nets;
    for(auto g : gates)
    {
        for(int n : g->nets)
        {
            nets.insert(n);
        }
    }
    for(auto p : pads)
    {
        for(int n : p->nets)
        {
            nets.insert(n);
        }
    }

    FILE* fid = fopen(filename, "w");
    if(fid != nullptr)
    {
        fprintf(fid, "%d %d\n", int(gates.size()), int(nets.size()));
        for(const auto& g : gates)
        {
            fprintf(fid, "%d %d", g->id, int(g->nets.size()));
            for(int ii : g->nets)
            {
                fprintf(fid, " %d", ii);
            }
            fprintf(fid, "\n");
        }
        fprintf(fid, "%d\n", int(pads.size()));
        for(auto p : pads)
        {
            fprintf(fid, "%d %d %d %d\n", p->id, p->nets[0], (int) p->x, (int) p->y);
        }
    }
}

std::vector<std::pair<int, int>> allPairs(const std::set<int>& in)
{
    std::vector<std::pair<int, int>> output;
    for(int ii : in)
    {
        for(int jj : in)
        {
            if(ii > jj)
            {
                output.emplace_back(ii, jj);
            }
        }
    }
    return output;
}

std::vector<Edge> reduce(const GateVector& gates, const GateVector& pads)
{
    // nets -> gates/pads
    std::map<int, std::set<int>> netsToGates;
    
    for(auto g : gates)
    {
        for(auto n : g->nets)
        {
            netsToGates[n].insert(g->id);
        }
    }
    for(auto p : pads)
    {
        for(auto n : p->nets)
        {
            netsToGates[n].insert(p->id);
        }
    }
    std::vector<Edge> output;
    for(const auto& s : netsToGates)
    {
        //int net = s.first;
        double weight = 2.0 / (s.second.size() * (s.second.size()-1));
        for(const auto& pair : allPairs(s.second))
        {
            Edge e;
            e.id1 = pair.first;
            e.id2 = pair.second;
            e.weight = weight;
            output.push_back(e);
        }
    }
    return output;
}

Placer::Placer(QString filename)
{
    readNetlist(filename.toLatin1().constData(), m_gates, m_pads);
    auto sr = std::make_shared<SubRegion>();
    sr->gates = m_gates;
    sr->pads = m_pads;
    m_subRegions.push_back(sr);
}

void Placer::placeGates(GateVector& gates, const GateVector& pads)
{
    int N = (int)gates.size();

    std::map<int, int> idToIdx;
    for(int ii = 0; ii < int(gates.size()); ii++)
    {
        const auto& g = gates[ii];
        idToIdx[g->id] = ii;
    }

    Eigen::VectorXd bx = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd by = Eigen::VectorXd::Zero(N);
    DMat A = DMat::Zero(N,N);

    std::vector<Edge> edges = reduce(gates,pads);
    for(const auto& e : edges)
    {
        bool isGate1 = (idToIdx.count(e.id1) > 0);
        bool isGate2 = (idToIdx.count(e.id2) > 0);
        if(isGate1 && isGate2)
        {
            // both are gates, add a connection
            A(idToIdx[e.id1], idToIdx[e.id2]) += e.weight;
            A(idToIdx[e.id2], idToIdx[e.id1]) += e.weight;
        }
        else if(!isGate1 && isGate2)
        {
            // First one is a pad, add to diagonal
            A(idToIdx[e.id2], idToIdx[e.id2]) += e.weight;
            auto pad = get(pads, e.id1);
            assert(pad->id == e.id1);
            bx(idToIdx[e.id2]) += e.weight * pad->x;
            by(idToIdx[e.id2]) += e.weight * pad->y;
        }
        else if(isGate1 && !isGate2)
        {
            // Second one is a pad, add to diagonal
            A(idToIdx[e.id1], idToIdx[e.id1]) += e.weight;
            auto pad = get(pads, e.id2);
            assert(pad->id == e.id2);
            bx(idToIdx[e.id1]) += e.weight * pad->x;
            by(idToIdx[e.id1]) += e.weight * pad->y;
        }
        else
        {
            // not possible, cannot directly connect two pads
            //assert(false);
        }
    }


    // The diagonal is the sum of each row
    for(int row = 0; row < N; row++)
    {
        double sum = 0.0;
        for(int col = 0; col < N; col++)
        {
            sum += A(row,col);
        }
        A(row,row) = sum;
    }

    // Flip the sign of everything except the diagonal
    for(int row = 0; row < N; row++)
    {
        for(int col = 0; col < N; col++)
        {
            if(row != col)
            {
                A(row, col) *= -1.0;
            }
        }
    }

    Eigen::VectorXd x(N);
    Eigen::VectorXd y(N);

    Eigen::ConjugateGradient<DMat, Eigen::Lower|Eigen::Upper> cg;
    cg.compute(A);
    x = cg.solve(bx);
    y = cg.solve(by);

    for(size_t ii = 0; ii < x.size(); ii++)
    {
        gates[ii]->x = x(ii);
        gates[ii]->y = y(ii);
    }
}

QVector<QPointF> Placer::getEdges(void)
{
    QVector<QPointF> output;
    auto edges = reduce(m_gates, m_pads);
    output.reserve(edges.size()*2);
    for(const auto& e : edges)
    {
        bool isGate1 = (1 <= e.id1) && (e.id1 <= m_gates.size());
        bool isGate2 = (1 <= e.id2) && (e.id2 <= m_gates.size());

        std::shared_ptr<Gate> g1;
        if(isGate1)
        {
            g1 = get(m_gates, e.id1);
        }
        else
        {
            g1 = get(m_pads, e.id1);
        }

        std::shared_ptr<Gate> g2;
        if(isGate2)
        {
            g2 = get(m_gates, e.id2);
        }
        else
        {
            g2 = get(m_pads, e.id2);
        }

        output.append(QPointF(g1->x/100.0,g1->y/100.0));
        output.append(QPointF(g2->x/100.0,g2->y/100.0));
    }
    return output;
}

void Placer::split(void)
{
    auto in = m_subRegions.front(); m_subRegions.pop_front();
    auto out1 = std::make_shared<SubRegion>();
    auto out2 = std::make_shared<SubRegion>();

    std::function<bool(std::shared_ptr<Gate>, std::shared_ptr<Gate>)> compare = nullptr;
    std::function<std::shared_ptr<Gate>(std::shared_ptr<Gate>)> clamp1 = nullptr;
    std::function<std::shared_ptr<Gate>(std::shared_ptr<Gate>)> clamp2 = nullptr;
    if(in->height() > in->width())
    {
        // split horizontal
        double midpoint = (in->top+in->bot) / 2.0;
        out1->bot   = in->bot;
        out1->top   = midpoint;
        out1->left  = in->left;
        out1->right = in->right;
        out2->bot   = midpoint;
        out2->top   = in->top;
        out2->left  = in->left;
        out2->right = in->right;

        compare = [](std::shared_ptr<Gate> a, std::shared_ptr<Gate> b)
        {
            if(a->y == b->y)
            {
                return a->x < b->x;
            }
            return a->y < b->y;
        };

        clamp1 = [midpoint](std::shared_ptr<Gate> a)
        {
            a->y = std::min(a->y, midpoint);
            return a;
        };

        clamp2 = [midpoint](std::shared_ptr<Gate> a)
        {
            a->y = std::max(a->y, midpoint);
            return a;
        };
    }
    else
    {
        // split vertical
        double midpoint = (in->left+in->right) / 2.0;
        out1->bot   = in->bot;
        out1->top   = in->top;
        out1->left  = in->left;
        out1->right = midpoint;
        out2->bot   = in->bot;
        out2->top   = in->top;
        out2->left  = midpoint;
        out2->right = in->right;

        compare = [](std::shared_ptr<Gate> a, std::shared_ptr<Gate> b)
        {
            if(a->x == b->x)
            {
                return a->y < b->y;
            }
            return a->x < b->x;
        };

        clamp1 = [midpoint](std::shared_ptr<Gate> a)
        {
            a->x = std::min(a->x, midpoint);
            return a;
        };

        clamp2 = [midpoint](std::shared_ptr<Gate> a)
        {
            a->x = std::max(a->x, midpoint);
            return a;
        };
    }

    // Sort the gates by position, smallest to largest.
    // Assign the lower half to out1 and the other half to out2.
    // In the case of an odd number, out1 gets an extra gate.

    std::sort(in->gates.begin(), in->gates.end(), compare);
    size_t half = (in->gates.size() + 1) / 2; // round up in case of odd number
    for(size_t ii = 0; ii < in->gates.size(); ii++)
    {
        if(ii < half)
        {
            out1->gates.push_back(clamp1(in->gates[ii]));
        }
        else
        {
            out2->gates.push_back(clamp2(in->gates[ii]));
        }
    }

    m_subRegions.push_back(out1);
    m_subRegions.push_back(out2);
}

void Placer::partition(void)
{
    auto sr = nextSubRegion();
    GateVector virtualPads = getVirtualPads(sr);
    placeGates(sr->gates, virtualPads);
    sr->placed = true;
    
    size_t maxGates = 0;
    for(auto s : m_subRegions)
    {
        maxGates = std::max(maxGates, s->gates.size());
    }
    qDebug() << "max gates per region =" << maxGates;
}

std::shared_ptr<SubRegion> Placer::nextSubRegion()
{
    // if no subregions ready to place, call split to make some more
    // return the next one
    if(m_subRegions.back()->placed)
    {
        split();
    }

    for(auto s : m_subRegions)
    {
        if(!s->placed)
        {
            return s;
        }
    }
    return nullptr;
}

GateVector Placer::getVirtualPads(std::shared_ptr<SubRegion> sr)
{
    GateVector output;
    output.reserve(m_gates.size() + m_pads.size());
    for(auto g : m_gates)
    {
        if(std::find(sr->gates.begin(), sr->gates.end(), g) == sr->gates.end())
        {
            auto vg = std::make_shared<Gate>(*g);
            vg->x = std::max(sr->left, std::min(sr->right, vg->x));
            vg->y = std::max(sr->bot,  std::min(sr->top,   vg->y));
            output.push_back(vg);
        }
    }
    for(auto p : m_pads)
    {
        auto vp = std::make_shared<Gate>(*p);
        vp->x = std::max(sr->left, std::min(sr->right, vp->x));
        vp->y = std::max(sr->bot,  std::min(sr->top,   vp->y));
        output.push_back(vp);
    }
    return output;
}

QVector<QPointF> Placer::getBorders(void)
{
    QVector<QPointF> output;
    for(auto s : m_subRegions)
    {
        output.append(QPointF(s->left/100.0, s->top/100.0));
        output.append(QPointF(s->right/100.0,s->top/100.0));

        output.append(QPointF(s->right/100.0,s->top/100.0));
        output.append(QPointF(s->right/100.0,s->bot/100.0));

        output.append(QPointF(s->left/100.0, s->bot/100.0));
        output.append(QPointF(s->right/100.0,s->bot/100.0));

        output.append(QPointF(s->left/100.0, s->bot/100.0));
        output.append(QPointF(s->left/100.0, s->top/100.0));
    }
    return output;
}
