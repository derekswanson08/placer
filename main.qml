import QtQuick 2.0
import QtQuick.Controls 2.12
import QtQuick.Window 2.12
import Placer 1.0


Item {
    id: m_main
    anchors.fill: parent

    Column
    {
        Row {
            height: 50
            id: m_menuBar
            Button {
                anchors.top: parent.top
                anchors.bottom: parent.bottom
                text: "Partition"
                onClicked: m_chipView.onPartition();
            }
            Button {
                anchors.top: parent.top
                anchors.bottom: parent.bottom
                text: "Reset"
                onClicked: m_chipView.onReset();
            }
        }

        Rectangle {
            height: 480
            width: 600
            Item {
                anchors.fill: parent
                MouseArea {
                    id: m_mouseArea
                    anchors.fill: parent
                    onWheel: {
                        if (wheel.modifiers & Qt.ControlModifier)
                        {
                            var zoomDelta = wheel.angleDelta.y / 120;
                            if(zoomDelta > 0)
                            {
                                // zooming in
                                m_chipView.zoomPosX += ((wheel.x / m_mouseArea.width)  - 0.5) / m_chipView.zoomLevel;
                                m_chipView.zoomPosY += (0.5 - (wheel.y / m_mouseArea.height)) / m_chipView.zoomLevel;
                            }
                            else
                            {
                                // zooming out. Don't change the zoom pos, unless we try to go
                                // past the minimum zoom level (1). In that case, reset the zoom pos
                                // to center.
                                if(m_chipView.zoomLevel == 1)
                                {
                                    m_chipView.zoomPosX = 0.5;
                                    m_chipView.zoomPosY = 0.5;
                                }
                            }

                            m_chipView.zoomLevel += zoomDelta;
                        }
                    }
                }

                ChipViewItem {
                    id: m_chipView
                    anchors.fill: parent
                    anchors.margins: 20
                }
            }
        }
    }
}
