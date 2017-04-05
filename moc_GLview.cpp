/****************************************************************************
** Meta object code from reading C++ file 'GLview.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "GLview.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'GLview.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_GLview_t {
    QByteArrayData data[18];
    char stringdata0[219];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_GLview_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_GLview_t qt_meta_stringdata_GLview = {
    {
QT_MOC_LITERAL(0, 0, 6), // "GLview"
QT_MOC_LITERAL(1, 7, 15), // "process_example"
QT_MOC_LITERAL(2, 23, 0), // ""
QT_MOC_LITERAL(3, 24, 7), // "inflate"
QT_MOC_LITERAL(4, 32, 11), // "randomNoise"
QT_MOC_LITERAL(5, 44, 6), // "smooth"
QT_MOC_LITERAL(6, 51, 7), // "sharpen"
QT_MOC_LITERAL(7, 59, 8), // "truncate"
QT_MOC_LITERAL(8, 68, 10), // "splitFaces"
QT_MOC_LITERAL(9, 79, 9), // "starFaces"
QT_MOC_LITERAL(10, 89, 14), // "splitLongEdges"
QT_MOC_LITERAL(11, 104, 18), // "collapseShortEdges"
QT_MOC_LITERAL(12, 123, 15), // "loopSubdivision"
QT_MOC_LITERAL(13, 139, 26), // "centerVerticesTangentially"
QT_MOC_LITERAL(14, 166, 9), // "flipEdges"
QT_MOC_LITERAL(15, 176, 4), // "crop"
QT_MOC_LITERAL(16, 181, 18), // "bilateralSmoothing"
QT_MOC_LITERAL(17, 200, 18) // "meshSimplification"

    },
    "GLview\0process_example\0\0inflate\0"
    "randomNoise\0smooth\0sharpen\0truncate\0"
    "splitFaces\0starFaces\0splitLongEdges\0"
    "collapseShortEdges\0loopSubdivision\0"
    "centerVerticesTangentially\0flipEdges\0"
    "crop\0bilateralSmoothing\0meshSimplification"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_GLview[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      16,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   94,    2, 0x08 /* Private */,
       3,    0,   95,    2, 0x08 /* Private */,
       4,    0,   96,    2, 0x08 /* Private */,
       5,    0,   97,    2, 0x08 /* Private */,
       6,    0,   98,    2, 0x08 /* Private */,
       7,    0,   99,    2, 0x08 /* Private */,
       8,    0,  100,    2, 0x08 /* Private */,
       9,    0,  101,    2, 0x08 /* Private */,
      10,    0,  102,    2, 0x08 /* Private */,
      11,    0,  103,    2, 0x08 /* Private */,
      12,    0,  104,    2, 0x08 /* Private */,
      13,    0,  105,    2, 0x08 /* Private */,
      14,    0,  106,    2, 0x08 /* Private */,
      15,    0,  107,    2, 0x08 /* Private */,
      16,    0,  108,    2, 0x08 /* Private */,
      17,    0,  109,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void GLview::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        GLview *_t = static_cast<GLview *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->process_example(); break;
        case 1: _t->inflate(); break;
        case 2: _t->randomNoise(); break;
        case 3: _t->smooth(); break;
        case 4: _t->sharpen(); break;
        case 5: _t->truncate(); break;
        case 6: _t->splitFaces(); break;
        case 7: _t->starFaces(); break;
        case 8: _t->splitLongEdges(); break;
        case 9: _t->collapseShortEdges(); break;
        case 10: _t->loopSubdivision(); break;
        case 11: _t->centerVerticesTangentially(); break;
        case 12: _t->flipEdges(); break;
        case 13: _t->crop(); break;
        case 14: _t->bilateralSmoothing(); break;
        case 15: _t->meshSimplification(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject GLview::staticMetaObject = {
    { &QOpenGLWidget::staticMetaObject, qt_meta_stringdata_GLview.data,
      qt_meta_data_GLview,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *GLview::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *GLview::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_GLview.stringdata0))
        return static_cast<void*>(const_cast< GLview*>(this));
    if (!strcmp(_clname, "QOpenGLFunctions"))
        return static_cast< QOpenGLFunctions*>(const_cast< GLview*>(this));
    return QOpenGLWidget::qt_metacast(_clname);
}

int GLview::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QOpenGLWidget::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 16)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 16;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 16)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 16;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
