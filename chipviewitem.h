#pragma once

#include <QtQuick/QQuickItem>
#include <memory>

class Placer;

class ChipViewItem : public QQuickItem
{
    Q_OBJECT

    Q_PROPERTY(QVector<QPointF> p READ points WRITE setPoints NOTIFY pointsChanged)
    Q_PROPERTY(QVector<QPointF> p READ borders WRITE setBorders NOTIFY bordersChanged)
    Q_PROPERTY(double zoomLevel READ zoomLevel WRITE setZoomLevel NOTIFY zoomLevelChanged)
    Q_PROPERTY(double zoomPosX READ zoomPosX WRITE setZoomPosX NOTIFY zoomPosXChanged)
    Q_PROPERTY(double zoomPosY READ zoomPosY WRITE setZoomPosY NOTIFY zoomPosYChanged)

public:
    ChipViewItem(QQuickItem *parent = 0);
    ~ChipViewItem();

    QSGNode *updatePaintNode(QSGNode *, UpdatePaintNodeData *);

    QVector<QPointF> points() const { return m_points; }
    void setPoints(const QVector<QPointF> p);

    QVector<QPointF> borders() const { return m_borders; }
    void setBorders(const QVector<QPointF> p);

    void setZoomLevel(double z);
    double zoomLevel(void) const { return m_zoomLevel; }

    void setZoomPosX(double x);
    double zoomPosX(void) const { return m_zoomPosX; }
    void setZoomPosY(double y);
    double zoomPosY(void) const { return m_zoomPosY; }

signals:
    void pointsChanged();
    void bordersChanged();
    void zoomLevelChanged();
    void zoomPosXChanged();
    void zoomPosYChanged();

public slots:
    void onPartition();
    void onReset();

private:
    QVector<QPointF> m_points;
    QVector<QPointF> m_borders;
    double m_zoomLevel = 1.0;
    double m_zoomPosX = 0.5;
    double m_zoomPosY = 0.5;
    std::shared_ptr<Placer> m_placer;
};
