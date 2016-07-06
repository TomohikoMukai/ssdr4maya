#ifndef RIGID_TRANSFORM_H
#define RIGID_TRANSFORM_H
#pragma once

#include <DirectXMath.h>

class RigidTransform
{
public:
    DirectX::XMFLOAT4A& Rotation()
    {
        return rotation;
    }
    const DirectX::XMFLOAT4A& Rotation() const
    {
        return rotation;
    }
    DirectX::XMFLOAT3A& Translation()
    {
        return translation;
    }
    const DirectX::XMFLOAT3A& Translation() const
    {
        return translation;
    }

public:
    RigidTransform()
        : rotation(0, 0, 0, 1.0f),
        translation(0, 0, 0)
    {
    }
    RigidTransform(const RigidTransform& src)
    {
        rotation = src.rotation;
        translation = src.translation;
    }
    RigidTransform(const DirectX::XMFLOAT4A& r, const DirectX::XMFLOAT3A& t)
    {
        rotation = r;
        translation = t;
    }
    explicit RigidTransform(const DirectX::XMFLOAT4X4A& m)
    {
        translation = DirectX::XMFLOAT3A(m._41, m._42, m._43);
        DirectX::XMMATRIX xm = DirectX::XMLoadFloat4x4A(&m);
        DirectX::XMVECTOR rv = DirectX::XMQuaternionRotationMatrix(xm);
        if (DirectX::XMVectorGetW(rv) < 0)
        {
            rv = DirectX::XMVectorNegate(rv);
        }
        DirectX::XMStoreFloat4A(&rotation, rv);
    }
    ~RigidTransform()
    {
    }

public:
    RigidTransform& operator =(const RigidTransform& src)
    {
        rotation = src.rotation;
        translation = src.translation;
        return *this;
    }
    RigidTransform& operator =(const DirectX::XMFLOAT4X4A& m)
    {
        translation = DirectX::XMFLOAT3A(m._41, m._42, m._43);
        DirectX::XMMATRIX xm = DirectX::XMLoadFloat4x4A(&m);
        DirectX::XMVECTOR rv = DirectX::XMQuaternionRotationMatrix(xm);
        if (DirectX::XMVectorGetW(rv) < 0)
        {
            rv = DirectX::XMVectorNegate(rv);
        }
        DirectX::XMStoreFloat4A(&rotation, rv);
        return *this;
    }
    void Set(const DirectX::XMFLOAT4A& r, const DirectX::XMFLOAT3A& t)
    {
        rotation = r;
        translation = t;
    }
    DirectX::XMVECTOR TransformCoord(DirectX::CXMVECTOR v) const
    {
        DirectX::XMVECTOR u = DirectX::XMVector3Rotate(v, DirectX::XMLoadFloat4A(&rotation));
        return DirectX::XMVectorAdd(DirectX::XMLoadFloat3A(&translation), u);
    }

public:
    static RigidTransform Identity()
    {
        return RigidTransform();
    }
    static RigidTransform Inverse(const RigidTransform& src)
    {
        RigidTransform is;
        DirectX::XMVECTOR xv = DirectX::XMLoadFloat4A(&src.rotation);
        DirectX::XMStoreFloat4A(&is.rotation, DirectX::XMQuaternionConjugate(xv));
        xv = DirectX::XMLoadFloat3A(&src.translation);
        DirectX::XMStoreFloat3A(&is.translation, DirectX::XMVectorNegate(xv));
        return is;
    }

public:
    DirectX::XMFLOAT4X4A ToMatrix4x4() const
    {
        DirectX::XMMATRIX m = DirectX::XMMatrixAffineTransformation(DirectX::XMVectorSet(1.0f, 1.0f, 1.0f, 1.0f),
            DirectX::XMVectorZero(),
            XMLoadFloat4A(&rotation),
            XMLoadFloat3A(&translation));
        DirectX::XMFLOAT4X4A r;
        DirectX::XMStoreFloat4x4A(&r, m);
        return r;
    }
    static RigidTransform FromMatrix4x4(const DirectX::XMFLOAT4X4A& m)
    {
        RigidTransform at;
        at.translation = DirectX::XMFLOAT3A(m._41, m._42, m._43);
        DirectX::XMMATRIX xm = DirectX::XMLoadFloat4x4A(&m);
        DirectX::XMVECTOR rv = DirectX::XMQuaternionRotationMatrix(xm);
        if (DirectX::XMVectorGetW(rv) < 0)
        {
            rv = DirectX::XMVectorNegate(rv);
        }
        DirectX::XMStoreFloat4A(&at.rotation, rv);
        return at;
    }

private:
    DirectX::XMFLOAT4A rotation;
    DirectX::XMFLOAT3A translation;
};

#endif //RIGID_TRANSFORM_H
