/****************************************************************************
** Meta object code from reading C++ file 'cmsc427.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.8.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../cmsc427.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'cmsc427.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.8.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_CMSC427Win_t {
    QByteArrayData data[5];
    char stringdata0[89];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_CMSC427Win_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_CMSC427Win_t qt_meta_stringdata_CMSC427Win = {
    {
QT_MOC_LITERAL(0, 0, 10), // "CMSC427Win"
QT_MOC_LITERAL(1, 11, 24), // "on_action_Quit_triggered"
QT_MOC_LITERAL(2, 36, 0), // ""
QT_MOC_LITERAL(3, 37, 24), // "on_action_Open_triggered"
QT_MOC_LITERAL(4, 62, 26) // "on_pushOpenOBJFile_clicked"

    },
    "CMSC427Win\0on_action_Quit_triggered\0"
    "\0on_action_Open_triggered\0"
    "on_pushOpenOBJFile_clicked"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_CMSC427Win[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   29,    2, 0x08 /* Private */,
       3,    0,   30,    2, 0x08 /* Private */,
       4,    1,   31,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void, QMetaType::Bool,    2,

       0        // eod
};

void CMSC427Win::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        CMSC427Win *_t = static_cast<CMSC427Win *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->on_action_Quit_triggered(); break;
        case 1: _t->on_action_Open_triggered(); break;
        case 2: _t->on_pushOpenOBJFile_clicked((*reinterpret_cast< bool(*)>(_a[1]))); break;
        default: ;
        }
    }
}

const QMetaObject CMSC427Win::staticMetaObject = {
    { &QMainWindow::staticMetaObject, qt_meta_stringdata_CMSC427Win.data,
      qt_meta_data_CMSC427Win,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *CMSC427Win::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *CMSC427Win::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_CMSC427Win.stringdata0))
        return static_cast<void*>(const_cast< CMSC427Win*>(this));
    return QMainWindow::qt_metacast(_clname);
}

int CMSC427Win::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QMainWindow::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
