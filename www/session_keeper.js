// Keeps the user's work if the page is closed by mistake. The files and settings are saved in this browser only
// (IndexedDB and localStorage on this computer); nothing is sent to the server. Closing the tab after loading a file
// asks for confirmation, and on the next visit a bar offers to restore the saved session or discard it.
(function () {
  const KEY = "electrofast_session_v1", DB = "electrofast_files", STORE = "files";
  let dirty = false, restoring = false;

  function openDB() {
    return new Promise((res, rej) => {
      const r = indexedDB.open(DB, 1);
      r.onupgradeneeded = () => r.result.createObjectStore(STORE);
      r.onsuccess = () => res(r.result); r.onerror = () => rej(r.error);
    });
  }
  function tx(db, mode, fn) {
    return new Promise((res, rej) => { const t = db.transaction(STORE, mode); const out = fn(t.objectStore(STORE));
      t.oncomplete = () => res(out && out.result !== undefined ? out.result : undefined); t.onerror = () => rej(t.error); });
  }
  async function putFiles(id, files) {
    const arr = [];
    for (const f of files) arr.push({ name: f.name, type: f.type, data: await f.arrayBuffer() });
    const db = await openDB(); await tx(db, "readwrite", s => s.put(arr, id));
  }
  async function getFiles(id) { const db = await openDB(); return tx(db, "readonly", s => s.get(id)); }
  async function clearAll() {
    try { localStorage.removeItem(KEY); } catch (e) {}
    try { const db = await openDB(); await tx(db, "readwrite", s => s.clear()); } catch (e) {}
  }
  function loadState() { try { return JSON.parse(localStorage.getItem(KEY)) || null; } catch (e) { return null; } }
  function saveState(st) { try { localStorage.setItem(KEY, JSON.stringify(st)); } catch (e) {} }

  // inputs that are events or table internals, not settings
  const skip = id => !id || id.startsWith(".clientdata") || /_cell_edit$|_rows_|_cells_|_search|_state$|_columns_selected$/.test(id);

  // settings: every change is saved
  $(document).on("shiny:inputchanged", function (ev) {
    if (restoring || skip(ev.name)) return;
    const el = document.getElementById(ev.name); if (!el) return;
    if (el.type === "file" || el.classList.contains("action-button") || el.tagName === "BUTTON") return;
    const st = loadState() || { inputs: {}, files: [] };
    st.inputs[ev.name] = ev.value; st.saved = new Date().toISOString(); saveState(st);
  });

  // files: saved when chosen
  $(document).on("change", "input[type=file]", async function () {
    if (restoring || !this.id || !this.files.length) return;
    dirty = true;
    try {
      await putFiles(this.id, this.files);
      const st = loadState() || { inputs: {}, files: [] };
      if (!st.files.includes(this.id)) st.files.push(this.id);
      st.saved = new Date().toISOString(); saveState(st);
    } catch (e) { console.warn("ElectroFast: could not save the file in this browser", e); }
  });

  // warning before the page is closed or reloaded with work in it
  window.addEventListener("beforeunload", function (e) { if (dirty) { e.preventDefault(); e.returnValue = ""; return ""; } });

  function waitForUpload(id, ms) {
    return new Promise(res => {
      const h = ev => { if (ev.name === id) { $(document).off("shiny:inputchanged", h); res(true); } };
      $(document).on("shiny:inputchanged", h); setTimeout(() => { $(document).off("shiny:inputchanged", h); res(false); }, ms);
    });
  }

  async function restore(st) {
    restoring = true;
    // files first (they set some choices automatically), then the saved settings on top
    for (const id of st.files || []) {
      const el = document.getElementById(id); const saved = await getFiles(id);
      if (!el || !saved || !saved.length) continue;
      const dt = new DataTransfer();
      for (const f of saved) dt.items.add(new File([f.data], f.name, { type: f.type }));
      el.files = dt.files; const done = waitForUpload(id, 60000);
      el.dispatchEvent(new Event("change", { bubbles: true })); await done;
    }
    for (const [id, val] of Object.entries(st.inputs || {})) {
      const el = document.getElementById(id); if (!el) continue;
      const b = $(el).data("shiny-input-binding"); if (!b || !b.receiveMessage) continue;
      try { b.receiveMessage(el, { value: val }); } catch (e) {}
    }
    restoring = false; dirty = (st.files || []).length > 0;
    // save the restored session again, so a second accidental close keeps it too; restored settings win over the
    // defaults recorded while the page was loading, and the files are still in this browser's database
    const cur = loadState() || { inputs: {}, files: [] };
    saveState({ inputs: Object.assign({}, cur.inputs, st.inputs), files: Array.from(new Set([...(cur.files || []), ...(st.files || [])])),
                saved: new Date().toISOString() });
  }

  function bar(st) {
    const when = st.saved ? new Date(st.saved).toLocaleString() : "earlier";
    const n = (st.files || []).length;
    const div = document.createElement("div");
    div.id = "electrofast-restore-bar";
    div.style.cssText = "position:fixed;left:0;right:0;bottom:0;z-index:2000;background:#1F2328;color:#fff;padding:10px 18px;" +
                        "display:flex;gap:12px;align-items:center;flex-wrap:wrap;font-size:0.92rem;";
    div.innerHTML = "<span>A previous session (" + n + " file input" + (n === 1 ? "" : "s") + ", saved " + when +
      ") is kept in this browser on this computer only.</span>" +
      "<button class='btn btn-sm btn-light' id='electrofast-restore'>Restore previous session</button>" +
      "<button class='btn btn-sm btn-outline-light' id='electrofast-discard'>Discard</button>";
    document.body.appendChild(div);
    // either choice removes the copy kept aside when the page opened
    const dropPending = () => { try { localStorage.removeItem(KEY + "_pending"); } catch (e) {} };
    document.getElementById("electrofast-restore").onclick = async () => { div.remove(); await restore(st); dropPending(); };
    document.getElementById("electrofast-discard").onclick = async () => { div.remove(); await clearAll(); dropPending(); };
  }

  $(document).one("shiny:connected", function () {
    const st = loadState();
    if (st && ((st.files && st.files.length) || Object.keys(st.inputs || {}).length)) {
      // keep the old copy aside until the user decides, so new changes do not overwrite it
      try { localStorage.setItem(KEY + "_pending", JSON.stringify(st)); localStorage.removeItem(KEY); } catch (e) {}
      if (st.files && st.files.length) bar(st); else { try { localStorage.removeItem(KEY + "_pending"); } catch (e) {} }
    }
  });
})();
