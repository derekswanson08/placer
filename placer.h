#pragma once

#include <QVector>
#include <QPointF>
#include <vector>
#include <memory>
#include <list>

struct SubRegion;
struct Gate;
typedef std::vector<std::shared_ptr<Gate>> GateVector;

class Placer
{
    public:
    Placer(QString filename);

    void partition();

    // Return the Rat's nest as a list of edges (point pairs)
    QVector<QPointF> getEdges(void);

    // Returns a list of edges (point pairs) that define the borders
    // between subregions
    QVector<QPointF> getBorders(void);

    private:
    void placeGates(GateVector& gates, const GateVector& pads);
    std::shared_ptr<SubRegion> nextSubRegion();
    GateVector getVirtualPads(std::shared_ptr<SubRegion> sr);
    void split();

    GateVector m_gates;
    GateVector m_pads;
    std::list<std::shared_ptr<SubRegion>> m_subRegions;
};

