#include <QtGui>
#include <QtQuick>

#include "chipviewitem.h"
#include "placer.h"

int main(int argc, char **argv)
{
    QGuiApplication app(argc, argv);

    qmlRegisterType<ChipViewItem>("Placer", 1, 0, "ChipViewItem");

    QQuickView view;
    QSurfaceFormat format = view.format();
    format.setSamples(16);
    view.setFormat(format);
    view.setSource(QUrl("qrc:///qml/main.qml"));
    view.show();

    return app.exec();
}
