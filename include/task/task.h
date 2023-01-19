#ifndef TASK_H_
#define TASK_H_

#include <cstdint>
#include <queue>
#include <vector>

enum TASK_TYPE {
    ADD_UINT,
    ADD_ULONG,
    ADD_DOUBLE,
    ADD_PAIR_DOUBLE_UINT,
    ADD_MIXED_PAIR_DOUBLE_UINT,
    ADD_UINT_WITH_REPLACE_PARENT,
    UINT_REPLACE_PARENT,
    MIN_UINT_WITH_PARENT,
    DIV_DOUBLE,
    DEC_UINT,
    DEC_ULONG,
    DEC_DOUBLE,
    DEC_PAIR_UINT_ULONG,
    SWAP_CIPHER_ENTRY
};

// When finished, the task result is written to the first operand.

#ifdef SGXBACKEND
#define DATA_SIZE 16
#define MAC_SIZE 16
#define NAC_SIZE 16

struct CipherEntry{
	uint8_t data[DATA_SIZE];	// Ciphertext content
	uint8_t nac[NAC_SIZE];		// This field store nonce + counter 
	uint8_t mac[MAC_SIZE];		// This field stores MAC of data entry fields

    // CipherEntry& operator=(const CipherEntry & ce) {
    //     memcpy(data, ce.data, DATA_SIZE);
    //     memcpy(nac, ce.nac, NAC_SIZE);
    //     memcpy(mac, ce.mac, MAC_SIZE);
    //     return *this;
    // }
};
typedef struct CipherEntry CipherEntry;
#endif

// #ifdef SSHEBACKEND
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>
#include <cstdlib>
#include <atomic>
#include <mutex>

class CryptoUtil;

#define MAX_CIPHER_TEXT_SIZE 850
#define MAX_PUBLIC_KEY_SIZE 1800

struct CipherEntry {
	uint8_t ct[MAX_CIPHER_TEXT_SIZE]; // HE ciphertext
	size_t tid; // tile index of the public key owner
    bool isShare;
    int64_t share[2];
    int plainNum;

    template <typename Archive> 
    void serialize(Archive &ar, const unsigned int version) 
    { 
        ar & ct;
        ar & tid; 
        ar & isShare;
        ar & share;
        ar & plainNum;
    }
};
typedef struct CipherEntry CipherEntry;

struct CipherEntryVec {
	std::vector<CipherEntry> ceVecs;

    template <typename Archive> 
    void serialize(Archive &ar, const unsigned int version) 
    { 
        ar & ceVecs;
    }
};
typedef struct CipherEntryVec CipherEntryVec;

struct PackedCipherEntry {
	uint8_t ct[MAX_CIPHER_TEXT_SIZE]; // HE ciphertext
	size_t tid; // tile index of the public key owner
    bool isShare;
    std::vector<int64_t> share;
    int plainNum;
    std::vector<int32_t> taskId;
    std::vector<int32_t> operandId;
    std::vector<int32_t> dstTid;

    template <typename Archive> 
    void serialize(Archive &ar, const unsigned int version) 
    { 
        ar & ct;
        ar & tid; 
        ar & isShare;
        ar & share;
        ar & plainNum;
        ar & taskId;
        ar & operandId;
        ar & dstTid;
    }
};
typedef struct PackedCipherEntry PackedCipherEntry;

struct ShareVec {
	std::vector<int64_t> shares;

    template <typename Archive> 
    void serialize(Archive &ar, const unsigned int version) 
    { 
        ar & shares;
    }
};
typedef struct ShareVec ShareVec;

struct PosVec {
    std::vector<uint64_t> pos;

    template <typename Archive> 
    void serialize(Archive &ar, const unsigned int version) 
    { 
        ar & pos;
    }
};
typedef struct PosVec PosVec;

struct ShareVecVec {
	std::vector<ShareVec> shareVecs;

    template <typename Archive> 
    void serialize(Archive &ar, const unsigned int version) 
    { 
        ar & shareVecs;
    }
};
typedef struct ShareVecVec ShareVecVec;

struct PublicKeyEntry {
    uint8_t pk[MAX_PUBLIC_KEY_SIZE];

    template <typename Archive> 
    void serialize(Archive &ar, const unsigned int version) 
    { 
        ar & pk;
    }
};
typedef struct PublicKeyEntry PublicKeyEntry;

// #endif

class TaskPayload {
public:
    bool isOperandEncrypted[2];
    virtual CipherEntry splitShareFromEncryptedOperand(int operandId0 = -1) = 0;
    virtual void writeShareToOperand(int64_t share, int operandId0 = -1, int operandId1 = -1) = 0;
    virtual ShareVec getOperandShare(int operandId0 = -1) = 0;
    virtual int getEncTid(int operandId0 = -1) = 0;
    virtual CipherEntry encryptShare(int operandId0 = -1, const uint32_t tid = -1) = 0;
    virtual void mergeEncryptedShare(CipherEntry& ce, int operandId0 = -1) = 0;
    virtual void writeEncryptedOperand(CipherEntry& ce, int operandId0 = -1) = 0;
    virtual CipherEntry* getCipherEntryPtr(int operandId0 = -1) = 0;
    virtual void unifyOperand() = 0;
    virtual void copyOperand(TaskPayload& srcTP, int srcOperandId0 = 1, int dstOperandId0 = 0) = 0;
    virtual void setCipherEntryPlainNum(int operandId0 = -1, int num = 0) = 0;
    static CryptoUtil& cryptoUtil;
    int operandNum;
    int plainNumPerCE;
    bool operandMask;
    bool useMask;
    
    TaskPayload () {
        isOperandEncrypted[0] = false;
        isOperandEncrypted[1] = false;
        operandMask = true;
        useMask = false;
    }
};

template<typename OperandType, typename IndexType>
class MinWithParent : public TaskPayload {
public:
    union {
        struct {
            OperandType plain_val;
            IndexType plain_parent;
        };
        CipherEntry enc;
    } operands[2];

    MinWithParent() {
        operandNum = 2;
        plainNumPerCE = 2;
    }

    CipherEntry splitShareFromEncryptedOperand(int operandId0) {
        std::vector<int64_t> shares;
        cryptoUtil.splitRandomShareFromCipherEntry(operands[operandId0].enc, 2, shares);
        for (int i=0; i<2; ++i) operands[operandId0].enc.share[i] = shares[i];
        // printf("split random share MWP %lld, %lld\n", shares[0], shares[1]);
        operands[operandId0].enc.isShare = true;
        if (operands[operandId0].enc.plainNum != 2) {
            printf("error MWP plain num %d\n", operands[operandId0].enc.plainNum);
            exit(-1);
        }
        return operands[operandId0].enc;
    }

    CipherEntry encryptShare(int operandId0, const uint32_t tid) {
        std::vector<int64_t> shares = {operands[operandId0].enc.share[0], operands[operandId0].enc.share[1]};
        CipherEntry ce = cryptoUtil.encryptToCipherEntry(shares, tid);
        return ce;
    }

    void mergeEncryptedShare(CipherEntry& ce, int operandId0) {
        std::vector<int64_t> shares = {operands[operandId0].enc.share[0], operands[operandId0].enc.share[1]};
        cryptoUtil.mergeShareIntoCipherEntry(ce, shares, ce.tid);
        // std::copy((uint8_t*)ce.ct, (uint8_t*)(ce.ct) + MAX_CIPHER_TEXT_SIZE, (uint8_t*)(operands[operandId0].enc.ct));
        operands[operandId0].enc = ce;
        operands[operandId0].enc.isShare = false;
        isOperandEncrypted[operandId0] = true;
        // operands[operandId0].enc.tid = ce.tid;
        return;
    }

    void writeShareToOperand(int64_t share, int operandId0, int operandId1) {
        operands[operandId0].enc.share[operandId1] = share;
        operands[operandId0].enc.isShare = true;
        isOperandEncrypted[operandId0] = true;
        return;
    }

    ShareVec getOperandShare(int operandId0) {
        ShareVec sv;
        sv.shares.push_back(operands[operandId0].enc.share[0]);
        sv.shares.push_back(operands[operandId0].enc.share[1]);
        return sv;
    }

    int getEncTid(int operandId0) {
        return operands[operandId0].enc.tid;
    }

    void writeEncryptedOperand(CipherEntry& ce, int operandId0) {
        operands[operandId0].enc = ce;
        operands[operandId0].enc.isShare = false;
        isOperandEncrypted[operandId0] = true;
    }

    CipherEntry* getCipherEntryPtr(int operandId0) {
        return &operands[operandId0].enc;
    }

    void unifyOperand() {
        operands[1] = operands[0];
        isOperandEncrypted[1] = isOperandEncrypted[0];
    }

    void copyOperand(TaskPayload& srcTP, int srcOperandId0, int dstOperandId0) {
        operands[dstOperandId0] = static_cast<MinWithParent&>(srcTP).operands[srcOperandId0];
    }

    void setCipherEntryPlainNum(int operandId0, int num) {
        operands[operandId0].enc.plainNum = num;
    }
};

template<typename OperandType>
class Add : public TaskPayload {
public:
    union {
        OperandType plain;
        CipherEntry enc;
    } operands[2];

    Add() {
        operandNum = 2;
        plainNumPerCE = 1;
    }

    CipherEntry splitShareFromEncryptedOperand(int operandId0) {
        int64_t share = cryptoUtil.splitRandomShareFromCipherEntry(operands[operandId0].enc);
        operands[operandId0].enc.share[0] = share;
        operands[operandId0].enc.isShare = true;
        return operands[operandId0].enc;
    }

    CipherEntry encryptShare(int operandId0, const uint32_t tid) {
        CipherEntry ce = cryptoUtil.encryptToCipherEntry(operands[operandId0].enc.share[0], tid);
        return ce;
    }

    void mergeEncryptedShare(CipherEntry& ce, int operandId0) {
        cryptoUtil.mergeShareIntoCipherEntry(ce, operands[operandId0].enc.share[0], ce.tid);
        operands[operandId0].enc = ce;
        operands[operandId0].enc.isShare = false;
        isOperandEncrypted[operandId0] = true;
        return;
    }

    void writeShareToOperand(int64_t share, int operandId0, int operandId1) {
        operands[operandId0].enc.share[0] = share;
        operands[operandId0].enc.isShare = true;
        isOperandEncrypted[operandId0] = true;
        return;
    }

    ShareVec getOperandShare(int operandId0) {
        ShareVec sv;
        sv.shares.push_back(operands[operandId0].enc.share[0]);
        sv.shares.push_back(operands[operandId0].enc.share[1]);
        return sv;
    }

    int getEncTid(int operandId0) {
        return operands[operandId0].enc.tid;
    }

    void writeEncryptedOperand(CipherEntry& ce, int operandId0) {
        operands[operandId0].enc = ce;
        operands[operandId0].enc.isShare = false;
        isOperandEncrypted[operandId0] = true;
    }

    CipherEntry* getCipherEntryPtr(int operandId0) {
        return &operands[operandId0].enc;
    }

    void unifyOperand() {
        operands[1] = operands[0];
        isOperandEncrypted[1] = isOperandEncrypted[0];
    }

    void copyOperand(TaskPayload& srcTP, int srcOperandId0, int dstOperandId0) {
        operands[dstOperandId0] = static_cast<Add&>(srcTP).operands[srcOperandId0];
    }

    void setCipherEntryPlainNum(int operandId0, int num) {
        operands[operandId0].enc.plainNum = num;
    }
};

template<typename OperandType, typename IndexType>
class AddWithReplaceParent : public MinWithParent<OperandType, IndexType> {
};

template<typename OperandType, typename IndexType>
class ReplaceParent : public MinWithParent<OperandType, IndexType> {
};

template<typename OperandType1, typename OperandType2>
class AddPair : public TaskPayload {
public:
    union {
        struct {
            OperandType1 plain_a;
            OperandType2 plain_b;
        };
        CipherEntry enc;
    } operands[2];

    AddPair() {
        operandNum = 2;
        plainNumPerCE = 2;
    }

    CipherEntry splitShareFromEncryptedOperand(int operandId0) {
        std::vector<int64_t> shares;
        cryptoUtil.splitRandomShareFromCipherEntry(operands[operandId0].enc, 2, shares);
        for (int i=0; i<2; ++i) operands[operandId0].enc.share[i] = shares[i];
        operands[operandId0].enc.isShare = true;
        return operands[operandId0].enc;
    }

    CipherEntry encryptShare(int operandId0, const uint32_t tid) {
        std::vector<int64_t> shares = {operands[operandId0].enc.share[0], operands[operandId0].enc.share[1]};
        CipherEntry ce = cryptoUtil.encryptToCipherEntry(shares, tid);
        return ce;
    }

    void mergeEncryptedShare(CipherEntry& ce, int operandId0) {
        std::vector<int64_t> shares = {operands[operandId0].enc.share[0], operands[operandId0].enc.share[1]};
        cryptoUtil.mergeShareIntoCipherEntry(ce, shares, ce.tid);
        operands[operandId0].enc = ce;
        operands[operandId0].enc.isShare = false;
        isOperandEncrypted[operandId0] = true;
        return;
    }

    void writeShareToOperand(int64_t share, int operandId0, int operandId1) {
        operands[operandId0].enc.share[operandId1] = share;
        operands[operandId0].enc.isShare = true;
        return;
    }

    ShareVec getOperandShare(int operandId0) {
        ShareVec sv;
        sv.shares.push_back(operands[operandId0].enc.share[0]);
        sv.shares.push_back(operands[operandId0].enc.share[1]);
        return sv;
    }

    int getEncTid(int operandId0) {
        return operands[operandId0].enc.tid;
    }

    void writeEncryptedOperand(CipherEntry& ce, int operandId0) {
        operands[operandId0].enc = ce;
        operands[operandId0].enc.isShare = false;
        isOperandEncrypted[operandId0] = true;
    }

    CipherEntry* getCipherEntryPtr(int operandId0) {
        return &operands[operandId0].enc;
    }

    void unifyOperand() {
        operands[1] = operands[0];
        isOperandEncrypted[1] = isOperandEncrypted[0];
    }

    void copyOperand(TaskPayload& srcTP, int srcOperandId0, int dstOperandId0) {
        operands[dstOperandId0] = static_cast<AddPair&>(srcTP).operands[srcOperandId0];
    }

    void setCipherEntryPlainNum(int operandId0, int num) {
        operands[operandId0].enc.plainNum = num;
    }
};

// The second member in the pair is plain
template<typename OperandType1, typename OperandType2>
class AddMixedPair : public TaskPayload {
public:
    struct {
        union {
            OperandType1 plain_a;
            CipherEntry enc_a;
        };
        OperandType2 plain_b;
    } operands[2];

    AddMixedPair() {
        operandNum = 2;
        plainNumPerCE = 1;
    }

    CipherEntry splitShareFromEncryptedOperand(int operandId0) {
        int64_t share = cryptoUtil.splitRandomShareFromCipherEntry(operands[operandId0].enc_a);
        operands[operandId0].enc_a.share[0] = share;
        operands[operandId0].enc_a.isShare = true;
        return operands[operandId0].enc_a;
    }

    CipherEntry encryptShare(int operandId0, const uint32_t tid) {
        CipherEntry ce = cryptoUtil.encryptToCipherEntry(operands[operandId0].enc_a.share[0], tid);
        return ce;
    }

    void mergeEncryptedShare(CipherEntry& ce, int operandId0) {
        cryptoUtil.mergeShareIntoCipherEntry(ce, operands[operandId0].enc_a.share[0], ce.tid);
        operands[operandId0].enc_a = ce;
        operands[operandId0].enc_a.isShare = false;
        isOperandEncrypted[operandId0] = true;
        return;
    }

    void writeShareToOperand(int64_t share, int operandId0, int operandId1) {
        operands[operandId0].enc_a.share[0] = share;
        operands[operandId0].enc_a.isShare = true;
        return;
    }

    ShareVec getOperandShare(int operandId0) {
        ShareVec sv;
        sv.shares.push_back(operands[operandId0].enc_a.share[0]);
        return sv;
    }

    int getEncTid(int operandId0) {
        return operands[operandId0].enc_a.tid;
    }

    void writeEncryptedOperand(CipherEntry& ce, int operandId0) {
        operands[operandId0].enc_a = ce;
        operands[operandId0].enc_a.isShare = false;
        isOperandEncrypted[operandId0] = true;
    }

    CipherEntry* getCipherEntryPtr(int operandId0) {
        return &operands[operandId0].enc_a;
    }

    void unifyOperand() {
        operands[1] = operands[0];
        isOperandEncrypted[1] = isOperandEncrypted[0];
    }

    void copyOperand(TaskPayload& srcTP, int srcOperandId0, int dstOperandId0) {
        operands[dstOperandId0] = static_cast<AddMixedPair&>(srcTP).operands[srcOperandId0];
    }

    void setCipherEntryPlainNum(int operandId0, int num) {
        operands[operandId0].enc_a.plainNum = num;
    }
};

template<typename OperandType>
class Div : public TaskPayload {
public:
    union {
        OperandType plain;
        CipherEntry enc;
    } operands[2];

    Div() {
        operandNum = 2;
        plainNumPerCE = 1;
    }

    CipherEntry splitShareFromEncryptedOperand(int operandId0) {
        int64_t share = cryptoUtil.splitRandomShareFromCipherEntry(operands[operandId0].enc);
        operands[operandId0].enc.share[0] = share;
        operands[operandId0].enc.isShare = true;
        return operands[operandId0].enc;
    }

    CipherEntry encryptShare(int operandId0, const uint32_t tid) {
        CipherEntry ce = cryptoUtil.encryptToCipherEntry(operands[operandId0].enc.share[0], tid);
        return ce;
    }

    void mergeEncryptedShare(CipherEntry& ce, int operandId0) {
        cryptoUtil.mergeShareIntoCipherEntry(ce, operands[operandId0].enc.share[0], ce.tid);
        operands[operandId0].enc = ce;
        operands[operandId0].enc.isShare = false;
        isOperandEncrypted[operandId0] = true;
        return;
    }

    void writeShareToOperand(int64_t share, int operandId0, int operandId1) {
        operands[operandId0].enc.share[0] = share;
        operands[operandId0].enc.isShare = true;
        return;
    }

    ShareVec getOperandShare(int operandId0) {
        ShareVec sv;
        sv.shares.push_back(operands[operandId0].enc.share[0]);
        return sv;
    }

    int getEncTid(int operandId0) {
        return operands[operandId0].enc.tid;
    }

    void writeEncryptedOperand(CipherEntry& ce, int operandId0) {
        operands[operandId0].enc = ce;
        operands[operandId0].enc.isShare = false;
        isOperandEncrypted[operandId0] = true;
    }

    CipherEntry* getCipherEntryPtr(int operandId0) {
        return &operands[operandId0].enc;
    }

    void unifyOperand() {
        operands[1] = operands[0];
        isOperandEncrypted[1] = isOperandEncrypted[0];
    }

    void copyOperand(TaskPayload& srcTP, int srcOperandId0, int dstOperandId0) {
        operands[dstOperandId0] = static_cast<Div&>(srcTP).operands[srcOperandId0];
    }

    void setCipherEntryPlainNum(int operandId0, int num) {
        operands[operandId0].enc.plainNum = num;
    }
};

class SwapCipherEntry : public TaskPayload {
public:
    union {
        CipherEntry enc;
    } operands[2];

    SwapCipherEntry() {
        operandNum = 2;
        plainNumPerCE = 2;
    }

    CipherEntry splitShareFromEncryptedOperand(int operandId0);

    CipherEntry encryptShare(int operandId0, const uint32_t tid);

    void mergeEncryptedShare(CipherEntry& ce, int operandId0);

    void writeShareToOperand(int64_t share, int operandId0, int operandId1);

    ShareVec getOperandShare(int operandId0);

    int getEncTid(int operandId0);

    void writeEncryptedOperand(CipherEntry& ce, int operandId0);

    CipherEntry* getCipherEntryPtr(int operandId0);

    void unifyOperand();

    void copyOperand(TaskPayload& srcTP, int srcOperandId0, int dstOperandId0);

    void setCipherEntryPlainNum(int operandId0, int num);
};

template<typename OperandType>
class DecryptEntry : public TaskPayload {
public:
    union {
        OperandType plain;
        CipherEntry enc;
    } operand;

    DecryptEntry() {
        operandNum = 1;
        plainNumPerCE = 1;
    }

    CipherEntry splitShareFromEncryptedOperand(int operandId0) {
        int64_t share = cryptoUtil.splitRandomShareFromCipherEntry(operand.enc);
        operand.enc.share[0] = share;
        operand.enc.isShare = true;
        return operand.enc;
    }

    CipherEntry encryptShare(int operandId0, const uint32_t tid) {
        CipherEntry ce = operand.enc;
        ce.tid = tid;
        // No encryption here
        return ce;
    }

    void mergeEncryptedShare(CipherEntry& ce, int operandId0) {
        operand.plain = cryptoUtil.template mergeShareAs<OperandType>(ce.share[0], operand.enc.share[0]);
        isOperandEncrypted[0] = false;
        return;
    }

    void writeShareToOperand(int64_t share, int operandId0, int operandId1) {
        operand.enc.share[0] = share;
        operand.enc.isShare = true;
        return;
    }

    ShareVec getOperandShare(int operandId0) {
        ShareVec sv;
        sv.shares.push_back(operand.enc.share[0]);
        return sv;
    }

    int getEncTid(int operandId0) {
        return operand.enc.tid;
    }

    void writeEncryptedOperand(CipherEntry& ce, int operandId0) {
        operand.enc = ce;
        operand.enc.isShare = false;
        isOperandEncrypted[0] = true;
    }

    CipherEntry* getCipherEntryPtr(int operandId0) {
        return &operand.enc;
    }

    void unifyOperand() {
        // Nothing to do
    }

    void copyOperand(TaskPayload& srcTP, int srcOperandId0, int dstOperandId0) {
        // Nothing to do
    }

    void setCipherEntryPlainNum(int operandId0, int num) {
        operand.enc.plainNum = num;
    }
};

template<typename OperandType1, typename OperandType2>
class DecryptEntryPair : public TaskPayload {
public:
    union {
        struct {
            OperandType1 plain_a;
            OperandType2 plain_b;
        };
        CipherEntry enc;
    } operand;

    DecryptEntryPair() {
        operandNum = 1;
        plainNumPerCE = 2;
    }

    CipherEntry splitShareFromEncryptedOperand(int operandId0) {
        std::vector<int64_t> shares;
        cryptoUtil.splitRandomShareFromCipherEntry(operand.enc, 2, shares);
        for (int i=0; i<2; ++i) operand.enc.share[i] = shares[i];
        operand.enc.isShare = true;
        return operand.enc;
    }

    CipherEntry encryptShare(int operandId0, const uint32_t tid) {
        CipherEntry ce = operand.enc;
        ce.tid = tid;
        // No encryption here
        return ce;
    }

    void mergeEncryptedShare(CipherEntry& ce, int operandId0) {
        std::vector<int64_t> shareVec0 = {ce.share[0], ce.share[1]};
        std::vector<int64_t> shareVec1 = {operand.enc.share[0], operand.enc.share[1]};
        cryptoUtil.mergeShareVec(shareVec0, shareVec1);
        operand.plain_a = cryptoUtil.template decodeFixedPointAs<OperandType1>(shareVec0[0]);
        operand.plain_b = cryptoUtil.template decodeFixedPointAs<OperandType1>(shareVec0[1]);
        isOperandEncrypted[operandId0] = false;
        return;
    }

    void writeShareToOperand(int64_t share, int operandId0, int operandId1) {
        operand.enc.share[operandId1] = share;
        operand.enc.isShare = true;
        return;
    }

    ShareVec getOperandShare(int operandId0) {
        ShareVec sv;
        sv.shares.push_back(operand.enc.share[0]);
        sv.shares.push_back(operand.enc.share[1]);
        return sv;
    }

    int getEncTid(int operandId0) {
        return operand.enc.tid;
    }

    void writeEncryptedOperand(CipherEntry& ce, int operandId0) {
        operand.enc = ce;
        operand.enc.isShare = false;
        isOperandEncrypted[0] = true;
    }

    CipherEntry* getCipherEntryPtr(int operandId0) {
        return &operand.enc;
    }

    void unifyOperand() {
        // Nothing to do
    }

    void copyOperand(TaskPayload& srcTP, int srcOperandId0, int dstOperandId0) {
        // Nothing to do
    }

    void setCipherEntryPlainNum(int operandId0, int num) {
        operand.enc.plainNum = num;
    }
};

struct Task {
    TASK_TYPE type;
    uint64_t vertexIndex; // Task destination vertex index
    bool finished;
    void* buf; // Task payload (content)

    uint64_t srcIndex = (uint64_t)-1; // Optional, task source vertex index

    uint64_t dstTid = 0; // Optional, task destination/cooperation tile index
    // For inter-graph edge, the task destination tile index is the owner of the destination vertex.
    // For intra-graph edge, the task destination is the encryptor of the source vertex data.
    uint64_t srcTid = 0; // Optional, task source tile index

    std::atomic<uint8_t> readyOperandNum{0}; // Optional

    bool isFinal = false;

    std::mutex mtx;

    bool isDummy = false;

    Task(const Task& t) {
        type = t.type;
        vertexIndex = t.vertexIndex;
        finished = t.finished;
        buf = t.buf;
        srcIndex = t.srcIndex;
        dstTid = t.dstTid;
        srcTid = t.srcTid;
        readyOperandNum = t.readyOperandNum.load();
        isFinal = t.isFinal;
        isDummy = t.isDummy;
    }

    Task& operator=(const Task& t) {
        type = t.type;
        vertexIndex = t.vertexIndex;
        finished = t.finished;
        buf = t.buf;
        srcIndex = t.srcIndex;
        dstTid = t.dstTid;
        srcTid = t.srcTid;
        readyOperandNum = t.readyOperandNum.load();
        isFinal = t.isFinal;
        isDummy = t.isDummy;
        return *this;        
    }

    Task() {
        type = (TASK_TYPE)0;
        vertexIndex = 0;
        finished = false;
        buf = NULL;
        srcIndex = 0;
        dstTid = 0;
        srcTid = 0;
        readyOperandNum = 0;
        isFinal = false;
        isDummy = false;
    }

    Task(TASK_TYPE a, uint64_t b, bool c, void* d, uint64_t e, uint64_t f, uint64_t g):
        type(a), vertexIndex(b), finished(c), buf(d), srcIndex(e), dstTid(f), srcTid(g), readyOperandNum(0), isFinal(false), isDummy(false) {
    }

    Task(TASK_TYPE a, uint64_t b, bool c, void* d, uint64_t e, uint64_t f, uint64_t g, uint64_t h, bool i):
        type(a), vertexIndex(b), finished(c), buf(d), srcIndex(e), dstTid(f), srcTid(g), readyOperandNum(h), isFinal(i), isDummy(false) {
    }

    Task(TASK_TYPE a, uint64_t b, bool c, void* d, uint64_t e, uint64_t f, uint64_t g, uint64_t h, bool i, bool isdm):
        type(a), vertexIndex(b), finished(c), buf(d), srcIndex(e), dstTid(f), srcTid(g), readyOperandNum(h), isFinal(i), isDummy(isdm) {
    }

#ifdef SSHEBACKEND
    template <typename Archive> 
    void serialize(Archive &ar, const unsigned int version) 
    { 
        ar & type;
        ar & vertexIndex; 
        ar & finished;
        ar & srcIndex;
        ar & dstTid;
        ar & srcTid;
        // ar & readyOperandNum;
        ar & isFinal;
    }
#endif

    void delete_task_content_buf() {
        switch (type) {
            case ADD_PAIR_DOUBLE_UINT:
                delete (AddPair<double, uint32_t>*)buf;
                break;
            case ADD_MIXED_PAIR_DOUBLE_UINT:
                delete (AddMixedPair<double, uint32_t>*)buf;
                break;
            case MIN_UINT_WITH_PARENT:
                delete (MinWithParent<uint32_t, uint64_t>*)buf;
                break;
            case ADD_DOUBLE:
                delete (Add<double>*)buf;
                break;
            case ADD_ULONG:
                delete (Add<uint64_t>*)buf;
                break;
            case ADD_UINT:
                delete (Add<uint32_t>*)buf;  
                break; 
            case ADD_UINT_WITH_REPLACE_PARENT:
                delete (AddWithReplaceParent<uint32_t, uint64_t>*)buf;  
                break;
            case UINT_REPLACE_PARENT:
                delete (ReplaceParent<uint32_t, uint64_t>*)buf;  
                break;     
            case DIV_DOUBLE:
                delete (Div<double>*)buf;
                break;
            case DEC_UINT:
                delete (DecryptEntry<uint32_t>*)buf;
                break;
            case DEC_ULONG:
                delete (DecryptEntry<uint64_t>*)buf;
                break;
            case DEC_DOUBLE:
                delete (DecryptEntry<double>*)buf;
                break;
            case DEC_PAIR_UINT_ULONG:
                delete (DecryptEntryPair<uint32_t, uint64_t>*)buf;
                break;
            case SWAP_CIPHER_ENTRY:
                delete (SwapCipherEntry*)buf;
                break;
            default:
                break;
        }
    }

    ~Task() {}
};

#endif