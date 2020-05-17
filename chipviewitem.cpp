#include <QtQuick/qsgnode.h>
#include <QtQuick/qsgvertexcolormaterial.h>
#include <algorithm>

#include "chipviewitem.h"
#include "placer.h"


ChipViewItem::ChipViewItem(QQuickItem *parent)
    : QQuickItem(parent)
{
    setClip(true);
    setFlag(ItemHasContents, true);
    onReset();
}

ChipViewItem::~ChipViewItem()
{
}

void ChipViewItem::onReset(void)
{
    m_placer = std::make_shared<Placer>(":/examples/fract");
    onPartition();
    setPoints(m_placer->getEdges());
    setBorders(m_placer->getBorders());
    setZoomLevel(1.0);
    setZoomPosX(0.5);
    setZoomPosY(0.5);
}

void ChipViewItem::onPartition(void)
{
    m_placer->partition();
    setPoints(m_placer->getEdges());
    setBorders(m_placer->getBorders());
}

void ChipViewItem::setPoints(const QVector<QPointF> p)
{
    if (p == m_points)
    {
        return;
    }

    m_points = p;
    emit pointsChanged();
    update();
}

void ChipViewItem::setBorders(const QVector<QPointF> p)
{
    if (p == m_borders )
    {
        return;
    }

    m_borders  = p;
    emit bordersChanged();
    update();
}

void ChipViewItem::setZoomLevel(double z)
{
    z = std::min(std::max(z, 1.0), 100.0);
    if(z == m_zoomLevel)
    {
        return;
    }

    m_zoomLevel = z;
    emit zoomLevelChanged();
    update();
}

void ChipViewItem::setZoomPosX(double x)
{
    x = std::min(std::max(x, 0.0), 1.0);
    if(x == m_zoomPosX)
    {
        return;
    }
    m_zoomPosX = x;
    emit zoomPosXChanged();
    update();
}

void ChipViewItem::setZoomPosY(double y)
{
    y = std::min(std::max(y, 0.0), 1.0);
    if(y == m_zoomPosY)
    {
        return;
    }
    m_zoomPosY = y;
    emit zoomPosYChanged();
    update();
}

QSGNode *ChipViewItem::updatePaintNode(QSGNode *oldNode, UpdatePaintNodeData *)
{
    QSGGeometryNode *node = nullptr;
    QSGGeometry *geometry = nullptr;
    QSGVertexColorMaterial *material = nullptr;

    int numVertices = m_points.size() + m_borders.size();
    if (!oldNode)
    {
        node = new QSGGeometryNode;
        geometry = new QSGGeometry(QSGGeometry::defaultAttributes_ColoredPoint2D(), 0);
        material = new QSGVertexColorMaterial;
    }
    else
    {
        node = static_cast<QSGGeometryNode *>(oldNode);
        geometry = node->geometry();
        material = static_cast<QSGVertexColorMaterial *>(node->material());
    }

    geometry->allocate(numVertices);
    geometry->setLineWidth(1);
    geometry->setDrawingMode(QSGGeometry::DrawLines);
    node->setGeometry(geometry);
    node->setFlag(QSGNode::OwnsGeometry);
    node->setMaterial(material);
    node->setFlag(QSGNode::OwnsMaterial);
    assert(geometry->vertexCount() == numVertices);

    QSizeF itemSize = size();
    QSGGeometry::ColoredPoint2D *vertices = geometry->vertexDataAsColoredPoint2D();
    for (int i = 0; i < m_points.size(); ++i) {
        QPointF pos = m_points.at(i);

        float x = (pos.x()-m_zoomPosX)*m_zoomLevel+.5;
        float y = (pos.y()-m_zoomPosY)*m_zoomLevel+.5;

        x = x * itemSize.width();
        y = (1.0-y) * itemSize.height();

        vertices[i].set(x, y, 0, 0, 255, 255);
    }
    for (int i = 0; i < m_borders.size(); ++i)
    {
        QPointF pos = m_borders.at(i);

        float x = (pos.x()-m_zoomPosX)*m_zoomLevel+.5;
        float y = (pos.y()-m_zoomPosY)*m_zoomLevel+.5;

        x = x * itemSize.width();
        y = (1.0-y) * itemSize.height();

        vertices[m_points.size()+i].set(x, y, 255, 0, 0, 255);
    }
    node->markDirty(QSGNode::DirtyGeometry);

    return node;
}

